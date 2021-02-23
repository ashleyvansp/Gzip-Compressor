/* uvgz.cpp

   Gzip-Compatible Lossless Data Compressor
   Ashley Van Spankeren

   Based on the starter code provided by B. Bird - 05/13/2020
*/

#include <iostream>
#include <vector>
#include <array>
#include <unordered_map>
#include <string>
#include "output_stream.hpp"
#include <algorithm>
#include <cassert>

// To compute CRC32 values, we can use this library
// from https://github.com/d-bahr/CRCpp
#define CRCPP_USE_CPP11
#include "CRC.h"
const int BLOCK_MAX = 100000;
const u32 EOB = 256;
std::vector<u32> type1_ll_lengths(288);
std::vector<u32> type1_ll_code(288);
std::vector<u32> type1_distance_lengths(32);
std::vector<u32> type1_distance_code(32);
std::vector<u32> base_length_offset_bits(6);
std::vector<u32> base_distances(14);
std::vector<u32> base_length_symbol_nums(28);
std::vector<u32> base_dist_symbol_nums(30);
const u32 window_size = 2500;

using Block = std::array<u8, BLOCK_MAX>;

//Given a vector of lengths where lengths.at(i) is the code length for symbol i,
//returns a vector V of unsigned int values, such that the lower lengths.at(i) bits of V.at(i)
//comprise the bit encoding for symbol i (using the encoding construction given in RFC 1951). Note that the encoding is in 
//MSB -> LSB order (that is, the first bit of the prefix code is bit number lengths.at(i) - 1 and the last bit is bit number 0).
//The codes for symbols with length zero are undefined.
std::vector< u32 > construct_canonical_code(std::vector<u32> const& lengths){

    unsigned int size = lengths.size();
    std::vector< unsigned int > length_counts(16,0); //Lengths must be less than 16 for DEFLATE
    u32 max_length = 0;
    for(auto i: lengths){
        //assert(i <= 15);
        length_counts.at(i)++;
        max_length = std::max(i, max_length);
    }
    length_counts[0] = 0; //Disregard any codes with alleged zero length

    std::vector< u32 > result_codes(size,0);

    //The algorithm below follows the pseudocode in RFC 1951
    std::vector< unsigned int > next_code(size,0);
    {
        //Step 1: Determine the first code for each length
        unsigned int code = 0;
        for(unsigned int i = 1; i <= max_length; i++){
            code = (code+length_counts.at(i-1))<<1;
            next_code.at(i) = code;
        }
    }
    {
        //Step 2: Assign the code for each symbol, with codes of the same length being
        //        consecutive and ordered lexicographically by the symbol to which they are assigned.
        for(unsigned int symbol = 0; symbol < size; symbol++){
            unsigned int length = lengths.at(symbol);
            if (length > 0)
                result_codes.at(symbol) = next_code.at(length)++;
        }  
    } 
    return result_codes;
}

// Given the length of the backreference, outputs the length code and offset (if applicable)
// For lengths less than or equal to 10 or equal to 258, there is no offset so we just get the code for that length and output it
// For any other length, it finds the base length and the offset and outputs them with the required number of bits
void write_ll_code_and_offset(u32 len, OutputBitStream& stream, const std::vector<u32>& ll_lengths, const std::vector<u32> ll_code){
    // Get the base length + offset amount
    u32 offset = 0;
    u32 offset_bits = 0;

    if (len <= 10){
        // No offset - just output the length
        // Length 3 = symbol 257, length 4 = symbol 258, etc.
        u32 symbol = len -  3 + 257;
        u32 length_code = ll_code.at(symbol);
        u32 length_code_bits = ll_lengths.at(symbol);
        // Output the base length code
        for(int idx = length_code_bits-1; idx >= 0; idx--)
            stream.push_bit((length_code>>(unsigned int)idx)&1);
    }
    else if (len == 258){
        // Also no offset - just output the length
        // I know this is hardcoded ...sry
        u32 symbol = 285;
        u32 length_code = ll_code.at(symbol);
        u32 length_code_bits = ll_lengths.at(symbol);
        // Output the base length code
        for(int idx = length_code_bits-1; idx >= 0; idx--)
            stream.push_bit((length_code>>(unsigned int)idx)&1);
    }
    else{
        while (len > base_length_offset_bits.at(offset_bits))
            offset_bits++;
        offset = (len - (base_length_offset_bits.at(offset_bits - 1) + 1)) % (1 << offset_bits);
        u32 base_len = len - offset;
        // Get symbol number
        u32 symbol = -1;
        for (int idx = 0; idx < 28; idx++){
            if (base_length_symbol_nums.at(idx) == base_len){
                symbol = idx + 257;
                break;
            }
        }
        if (symbol == -1){
            std::cout << "No matched symbol for base length " << base_len << "\n";
            std::abort();
        }

        u32 base_len_code = ll_code.at(symbol);
        u32 length_code_bits = ll_lengths.at(symbol);

        // Output the base length code
        for(int idx = length_code_bits-1; idx >= 0; idx--)
            stream.push_bit((base_len_code>>(unsigned int)idx)&1);
        // Output the offset - RULE #1
        stream.push_bits(offset, offset_bits);
    }
}

// Given the distance of the backreference, outputs the distance code and offset (if applicable)
void write_dist_code_and_offset(u32 dist, OutputBitStream& stream, const std::vector<u32>& dist_codes, const std::vector<u32>& dist_lengths){
    if (dist < 5){
        // No offset - just output the distance
        u32 code = dist_codes.at(dist - 1);
        u32 bits = dist_lengths.at(dist - 1);
        for(int idx = bits - 1; idx >= 0; idx--)
            stream.push_bit((code>>(unsigned int)idx) & 1);
    }else{
        u32 offset_bits = 1;
        while (dist > 1 << (offset_bits + 2))
            offset_bits++;
        u32 base_dist;
        u32 symbol = -1;
        if (dist > 24576){
            symbol = 29;
            base_dist = 24577;
        }else{
            for (u32 idx = 0; idx < 30; idx++){
                if (base_dist_symbol_nums.at(idx) > dist){
                    symbol = idx - 1;
                    base_dist = base_dist_symbol_nums.at(idx - 1);
                    break;
                }
            }
        }
        if (symbol == -1){
            std::cout << "No base distance found for distance " << dist << "\n";
            std::abort();
        }

        u32 offset = dist - base_dist;
        u32 base_dist_code = dist_codes.at(symbol);
        u32 bits = dist_lengths.at(symbol);
        // Output the base distance
        for(int idx = bits - 1; idx >= 0; idx--)
            stream.push_bit((base_dist_code>>(unsigned int)idx)&1);
        // Output the offset - RULE #1
        stream.push_bits(offset, offset_bits);
    }
}

void generate_type1(OutputBitStream& stream, Block& block_data, u32 block_size, bool is_last){
    stream.push_bit(is_last ? 1 : 0);
    stream.push_bits(1, 2);

    // Try with backreferences
    // First try: start furthest away. if a backreference is found, output it right away
    // upgrades: continue looking for backreferences to find one that's closer or longer

    for (u32 cur = 0; cur < block_size; cur++){
        u32 symbol = block_data.at(cur);
        u32 bits = type1_ll_lengths.at(symbol);
        u32 code = type1_ll_code.at(symbol);
        bool found_ref = false; 
        // Try sliding window of window size 
        u32 window = cur > window_size ? cur - window_size : 0;
        for (u32 hist = window; hist < cur; hist++){
            if (block_data.at(hist) == block_data.at(cur)){
                // Found a backreference; now find its length
                u32 len = 1;
                u32 dist = cur - hist;
                while (len < 258 && cur + len < block_size && block_data.at(hist + len) == block_data.at(cur + len)){
                    len++;
                }
                // Only care about backreferences longer than 3
                if (len > 2){
                    found_ref = true;
                    write_ll_code_and_offset(len, stream, type1_ll_lengths, type1_ll_code);
                    write_dist_code_and_offset(dist, stream, type1_distance_code, type1_distance_lengths);

                    // Now move on to look for backreferences at the next available spot
                    cur += (len - 1); // k will get incremented in the outer for loop
                    hist = cur + 1;
                }
            }
        }

        // If no backreference was found, output the literal
        if (! found_ref){
        for(int i = bits-1; i >= 0; i--)
            stream.push_bit((code>>(unsigned int)i)&1);
        }
    }

    // Append the EOB marker
    {
        u32 symbol = EOB;
        u32 bits = type1_ll_lengths.at(symbol);
        u32 code = type1_ll_code.at(symbol);
        for(int i = bits-1; i >= 0; i--)
            stream.push_bit((code>>(unsigned int)i)&1);
    }
}

// Used when finding frequencies in the initial stages of creating a type 2 block
// Gets called from type2_LZSS_frequencies when a backreference is found
// Finds the base distance and increments its frequency
// If the base distance hadn't been seen before, dist_used is incremented too (dist_used is helpful for eventually finding HDIST)
void get_dist_freq(u32 dist, std::vector<u32>& dist_counts, u32& dist_used){
    if (dist < 5){
        // No offset - just output the distance
        if (dist_counts.at(dist - 1) == 0)
            dist_used++;
        dist_counts.at(dist - 1)++;
    }else{
        u32 offset_bits = 1;
        while (dist > 1 << (offset_bits + 2))
            offset_bits++;
        u32 base_dist;
        u32 symbol = -1;
        if (dist > 24576){
            symbol = 29;
            base_dist = 24577;
        }else{
            for (u32 idx = 0; idx < 30; idx++){
                if (base_dist_symbol_nums.at(idx) > dist){
                    symbol = idx - 1;
                    base_dist = base_dist_symbol_nums.at(idx - 1);
                    break;
                }
            }
        }
        if (symbol == -1){
            std::cout << "No base distance found for distance " << dist << "\n";
            std::abort();
        }

        // Output the base distance
        if (dist_counts.at(symbol) == 0)
            dist_used++;
        dist_counts.at(symbol)++;
    }
}

// Used when finding frequencies in the initial stages of creating a type 2 block
// Gets called from type2_LZSS_frequencies when a backreference is found
// Finds the base length and increments its frequency
// If the base length hadn't been seen before, ll_used is incremented too (ll_used is helpful for eventually finding HLIT)
void get_ll_freq(u32 len, std::vector<u32>& ll_counts, u32& ll_used){
    // Get the base length + offset amount
    u32 offset = 0;
    u32 offset_bits = 0;

    if (len <= 10){
        // No offset - just output the length
        // Length 3 = symbol 257, length 4 = symbol 258, etc.
        u32 symbol = len -  3 + 257;

        // Update the base length frequency
        if (ll_counts.at(symbol) == 0)
            ll_used++;
        ll_counts.at(symbol)++;
    }
    else if (len == 258){
        // Also no offset - just output the length
        // I know this is hardcoded ...sry
        u32 symbol = 285;
        // Update the frequency of the base length
        if (ll_counts.at(symbol) == 0)
            ll_used++;
        ll_counts.at(symbol)++;
    }
    else{
        while (len > base_length_offset_bits.at(offset_bits))
            offset_bits++;
        offset = (len - (base_length_offset_bits.at(offset_bits - 1) + 1)) % (1 << offset_bits);
        u32 base_len = len - offset;
        // Get symbol number
        u32 symbol = -1;
        for (int idx = 0; idx < 28; idx++){
            if (base_length_symbol_nums.at(idx) == base_len){
                symbol = idx + 257;
                break;
            }
        }
        if (symbol == -1){
            std::cout << "No matched symbol for base length " << base_len << "\n";
            std::abort();
        }

        // Update the base length frequency
        if (ll_counts.at(symbol) == 0)
            ll_used++;
        ll_counts.at(symbol)++;
    }
}

// Runs LZSS on the block contents to determine the frequency of each length/ literal/ distance symbol
void type2_LZSS_frequencies(Block& block_data, u32 block_size, std::vector<u32>& ll_counts, std::vector<u32>& dist_counts, u32& ll_used, u32& dist_used){
    // Try with backreferences
    // First try: start furthest away. if a backreference is found, output it right away
    // upgrades: continue looking for backreferences to find one that's closer or longer

    for (u32 cur = 0; cur < block_size; cur++){
        u32 symbol = block_data.at(cur);
     //   u32 bits = type1_ll_lengths.at(symbol);
     //   u32 code = type1_ll_code.at(symbol);
        bool found_ref = false; 
        // Try sliding window of window size 
        u32 window = cur > window_size ? cur - window_size : 0;
        for (u32 hist = window; hist < cur; hist++){
            if (block_data.at(hist) == block_data.at(cur)){
                // Found a backreference; now find its length
                u32 len = 1;
                u32 dist = cur - hist;
                while (len < 258 && cur + len < block_size && block_data.at(hist + len) == block_data.at(cur + len)){
                    len++;
                }
                // Only care about backreferences longer than 3
                if (len > 2){
                    found_ref = true;
                    get_ll_freq(len, ll_counts, ll_used);
                    get_dist_freq(dist, dist_counts, dist_used);

                    // Now move on to look for backreferences at the next available spot
                    cur += (len - 1); // k will get incremented in the outer for loop
                    hist = cur + 1;
                }
            }
        }

        // If no backreference was found, output the literal
        if (! found_ref){
            ll_counts.at(symbol)++;
            if (ll_counts.at(symbol) == 1)
                ll_used++;
        }
    }

    // Append the EOB marker
    {
        u32 symbol = EOB;
        ll_counts.at(symbol)++;
        if (ll_counts.at(symbol) == 1)
            ll_used++;
    }
}

// Runs LZSS on the block_data and outputs directly to stream using codes from the ll_code and dist_code
void type2_LZSS(OutputBitStream& stream, Block& block_data, u32 block_size, const std::vector<u32>& ll_code, const std::vector<u32>& ll_lengths, const std::vector<u32>& dist_code, const std::vector<u32>& dist_lengths){

    // Try with backreferences
    // First try: start furthest away. if a backreference is found, output it right away
    // upgrades: continue looking for backreferences to find one that's closer or longer

    for (u32 cur = 0; cur < block_size; cur++){
        u32 symbol = block_data.at(cur);
        u32 bits = ll_lengths.at(symbol);
        u32 code = ll_code.at(symbol);
        bool found_ref = false; 
        // Try sliding window of window_size 
        u32 window = cur > window_size ? cur - window_size : 0;
        for (u32 hist = window; hist < cur; hist++){
            if (block_data.at(hist) == block_data.at(cur)){
                // Found a backreference; now find its length
                u32 len = 1;
                u32 dist = cur - hist;
                while (len < 258 && cur + len < block_size && block_data.at(hist + len) == block_data.at(cur + len)){
                    len++;
                }
                // Only care about backreferences longer than 3
                if (len > 2){
                    found_ref = true;
                    
                    write_ll_code_and_offset(len, stream, ll_lengths, ll_code);
                    write_dist_code_and_offset(dist, stream, dist_code, dist_lengths);
                    // Now move on to look for backreferences at the next available spot
                    cur += (len - 1); // k will get incremented in the outer for loop
                    hist = cur + 1;
                }
            }
        }

        // If no backreference was found, output the literal (or increment the literal's frequency)
        if (! found_ref){
            for(int i = bits-1; i >= 0; i--)
                stream.push_bit((code>>(unsigned int)i)&1);
        }
    }

    // Append the EOB marker
    {
        u32 symbol = EOB;
        u32 bits = ll_lengths.at(symbol);
        u32 code = ll_code.at(symbol);
        for(int i = bits-1; i >= 0; i--)
            stream.push_bit((code>>(unsigned int)i)&1);

    }
}

// Returns true if and only if the Kraft-McMillan inequality is met
bool KM_inequality(std::vector<u32> &W){
    float sum = 0;
    for (int i = 0; i < W.size(); i++){
        sum += (1.0 / (1 << W.at(i)));
    }
    return sum == 1.0;
}

// Given a vector containing the lengths of prefix codes, changes the lengths until the KM inequality is met
void make_optimal(std::vector<u32>& W){
    // Find the longest code length (bottom of tree)
    // Decrement its length
    // Does it meet the inequality?
    // Rinse repeat

    while (! KM_inequality(W)){
        u32 max_idx = std::max_element(W.begin(), W.end()) - W.begin();
        W.at(max_idx)--;
    }
}

// Based on pseudocode documented in:
// https://people.eng.unimelb.edu.au/ammoffat/abstracts/compsurv19moffat.pdf
// Assigns nonzero code lengths to zero frequencies
// Returns true if and only if the lengths equal the KM inequality
void calc_huff_lens(std::vector<u32>& W, int n){
    // Phase 1
    int leaf = n - 1;
    int root = n - 1;
    int next;
    for (next = n - 1; next > 0; next--){
        // Find first child
        if (leaf < 0 || (root > next && W.at(root) < W.at(leaf))){
            // Use internal node
            W.at(next) = W.at(root);
            W.at(root) = next;
            root--;
        }else{
            // Use leaf node
            W.at(next) = W.at(leaf);
            leaf--;
        }
        // Find second child
        // Repeat the above, but adding to W.at(next) rather than assigning to it
        if (leaf < 0 || (root > next && W.at(root) < W.at(leaf))){
            // Use internal node
            W.at(next) += W.at(root);
            W.at(root) = next;
            root--;
        }else{
            // Use leaf node
            W.at(next) += W.at(leaf);
            leaf--;
        }
    }

    // Phase 2
    W.at(1) = 0;
    for (next = 2; next < n; next++){
        W.at(next) = W.at(W.at(next)) + 1;
    }
    // Phase 3
    u32 avail = 1;
    u32 used = 0;
    u32 depth = 0;
    root = 1;
    next = 0;
    while (avail > 0){
        // Count internal nodes used at depth 'depth'
        while (root < n && W.at(root) == depth){
            used++;
            root++;
        }
        // Assign as leaves any nodes that are not internal
        while (avail > used){
            W.at(next) = depth;
            next++;
            avail--;
        }
        // Move to next depth
        avail = 2 * used;
        depth++;
        used = 0;
    }

    // W.at(i) now contains the length l_i of the ith codeword
    // Needs to meet the KM inequality!
    make_optimal(W);
}

// Takes as input a vector of code lengths
// Applies RLE using the 19-symbol code length encoding
// Returns a vector of the code length encoding
std::vector<u32> rle(std::vector<u32> code_lengths){
    std::vector<u32> cl_data(code_lengths.size());
    u16 len = code_lengths.size();
    u16 run_val = code_lengths.at(0);
    u16 run_length = 0;
    u16 cl_ind = 0;
    u16 i = 0;
    while (i < len){
        run_val = code_lengths.at(i);
        run_length = 0;
        // Matches no matter what for one iteration, so run length > 0 after the loop
        while (i < len && code_lengths.at(i) == run_val){
            run_length++;
            i++;
        }
        if (run_length < 3){
            // No notable run
            for (int k = 0; k < run_length; k++){
                cl_data.at(cl_ind++) = run_val;
            }
        } else if (run_val == 0 && run_length < 11){
            // Short run of zeroes
            cl_data.at(cl_ind++) = 17;
            cl_data.at(cl_ind++) = run_length - 3;
        } else if (run_val == 0 && run_length < 139){
            // Run of zeroes that can be encoded with one symbol
            cl_data.at(cl_ind++) = 18;
            cl_data.at(cl_ind++) = run_length - 11;
        } else if (run_val == 0){
            // Long run of zeroes; needs multiple symbols
            cl_data.at(cl_ind++) = 18;
            cl_data.at(cl_ind++) = 127;
            run_length-=138;
            // Leftover zeroes
            if (run_length < 3){
                for (int k = 0; k < run_length; k++){
                    cl_data.at(cl_ind++) = 0;
                }
            }else if (run_length < 11){
                cl_data.at(cl_ind++) = 17;
                cl_data.at(cl_ind++) = run_length - 3;
            }else{
                cl_data.at(cl_ind++) = 18;
                cl_data.at(cl_ind++) = run_length - 11;
            }
        } else{
            // Run of non zero values
            // Something is off with my 16s so for now just omit
            for (int k = 0; k < run_length; k++){
                cl_data.at(cl_ind++) = run_val;
            }
        }
    }
    cl_data.resize(cl_ind);
    return cl_data;
}

// Takes cl_symbols as the RLE of the code lengths (from above) and writes them to stream using the cl_codes and their lengths
void write_code_length_encoding(OutputBitStream& stream, const std::vector<u32>& cl_symbols, const std::vector<u32>& cl_codes, const std::vector<u32>& code_lengths){

    for (int k = 0; k < cl_symbols.size(); k++){
        u32 symbol = cl_symbols.at(k);
        u32 bits = code_lengths.at(symbol);
        u32 code = cl_codes.at(symbol);

        if (symbol == 16){
            // First write the CL code for symbol 16
            for(int i = bits-1; i >= 0; i--){
                stream.push_bit((code>>(unsigned int)i)&1);
            }
            // Then write the length with 2 bits - RULE #1
            stream.push_bits(cl_symbols.at(k + 1), 2);
            k++;
        }else if (symbol == 17){
            // First write the CL code for symbol 17
            for(int i = bits-1; i >= 0; i--){
                stream.push_bit((code>>(unsigned int)i)&1);
            }

            // Then write the length with 3 bits - RULE #1
            stream.push_bits(cl_symbols.at(k + 1), 3);
            k++;
        }else if (symbol == 18){
            // First write the CL code for symbol 18
            for(int i = bits-1; i >= 0; i--){
                stream.push_bit((code>>(unsigned int)i)&1);
            }
            // Then write the length with 7 bits - RULE #1
            stream.push_bits(cl_symbols.at(k + 1), 7);
            k++;
        }else{
            for(int i = bits-1; i >= 0; i--){
                stream.push_bit((code>>(unsigned int)i)&1);
            }
        }
    }
}

// Input: the RLE of the LL code lengths and distance code lengths
// Output: a vector V where V.at(i) is the number of times i occurs in the input
std::vector<u32> get_cl_freqs(std::vector<u32> ll_rle, std::vector<u32> dist_rle, u32 &cl_used){
    std::vector<u32> cl_freqs(19,0);
    cl_used = 0;
    for (int i = 0; i < ll_rle.size(); i++){
        u32 symbol = ll_rle.at(i);
        if (cl_freqs.at(symbol) == 0)
            cl_used++;
        if (symbol > 15)
            i++;
        cl_freqs.at(symbol)++;

    }
    for (int i = 0; i < dist_rle.size(); i++){
        u32 symbol = dist_rle.at(i);
        if (cl_freqs.at(symbol) == 0)
            cl_used++;
        if (symbol > 15)
            i++;

        cl_freqs.at(symbol)++;
    }
    return cl_freqs;
}

void write_cl_data(OutputBitStream& stream, std::vector<u32> const& ll_code_lengths, std::vector<u32> const& dist_code_lengths ){

    //Variables are named as in RFC 1951
    assert(ll_code_lengths.size() >= 257); //There needs to be at least one use of symbol 256, so the ll_code_lengths table must have at least 257 elements
    unsigned int HLIT = ll_code_lengths.size() - 257;

    unsigned int HDIST = 0;
    if (dist_code_lengths.size() == 0){
        //Even if no distance codes are used, we are required to encode at least one.
    }else{
        HDIST = dist_code_lengths.size() - 1;
    }

    std::vector<u32> ll_rle = rle(ll_code_lengths);
    std::vector<u32> dist_rle = rle(dist_code_lengths);
    u32 cl_used = 0;
    std::vector<u32> cl_freqs = get_cl_freqs(ll_rle, dist_rle, cl_used);
    std::vector<u32> cl_nz_freqs(cl_used, 0);
    std::vector<u32> cl_nz_symbols(cl_used, 0);
    // Get rid of zeroes in the cl frequency table so we can get the code lengths
    u32 cl_used_2 = 0;
    for (int i = 0; i < 19; i++){
        if (cl_freqs.at(i) > 0){
            cl_nz_freqs.at(cl_used_2) = cl_freqs.at(i);
            cl_nz_symbols.at(cl_used_2++) = i;
        }
    }
    // Get cl code lengths
    calc_huff_lens(cl_nz_freqs, cl_nz_freqs.size());

    // Put it back to original size
    std::vector<u32> cl_code_lengths(19, 0);
    for (int k = 0; k < cl_code_lengths.size(); k++){
        std::vector<u32>::iterator it = find (cl_nz_symbols.begin(), cl_nz_symbols.end(), k);
        if (it != cl_nz_symbols.end()){
            u32 index = distance(cl_nz_symbols.begin(), it);
            u32 code_length = cl_nz_freqs.at(index);
            cl_code_lengths.at(k) = code_length;
            if (code_length > 7){
                cl_code_lengths = {4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 5};
                break;
            }
        }else{
            cl_code_lengths.at(k) = 0;
        }
    }

    // Get the CL code
    std::vector<u32> cl_code = construct_canonical_code(cl_code_lengths);

    unsigned int HCLEN = 18;
    
    //The lengths are written in a strange order, dictated by RFC 1951
    //(This seems like a sadistic twist of the knife, but there is some amount of weird logic behind the ordering)
    std::vector<u32> cl_permutation {16, 17, 18, 0, 8, 7, 9, 6, 10, 5, 11, 4, 12, 3, 13, 2, 14, 1, 15};

    while (cl_code_lengths.at(cl_permutation.at(HCLEN)) == 0)
        HCLEN--;
    HCLEN-=3;

    //Push HLIT, HDIST and HCLEN. These are all numbers so Rule #1 applies
    stream.push_bits(HLIT, 5);
    stream.push_bits(HDIST,5);
    stream.push_bits(HCLEN,4);

    //Now push each CL code length in 3 bits (the lengths are numbers, so Rule #1 applies)
    //assert(HCLEN+4 == 19);
    for (unsigned int i = 0; i < HCLEN+4; i++)
        stream.push_bits(cl_code_lengths.at(cl_permutation.at(i)),3);

    write_code_length_encoding(stream, ll_rle, cl_code, cl_code_lengths);
    write_code_length_encoding(stream, dist_rle, cl_code, cl_code_lengths);

}

// First attempt - brute force
// 1. Run "fake" LZSS (ie. no output) to get literal / length / distance counts
// 2. Construct prefix codes based on frequencies of symbols
// 3. Output header now that we have all the necessary info
// 4. Run LZSS again but with output this time
// Returns true if there was success
// Returns false if failure such as code length greater than 15
bool generate_type2(OutputBitStream& stream, Block& block_data, u32 block_size, bool is_last){


    std::vector<u32> ll_counts(288, 0);
    std::vector<u32> dist_counts(32, 0);
    u32 ll_used = 0;
    u32 dist_used = 0;
    type2_LZSS_frequencies(block_data, block_size, ll_counts, dist_counts, ll_used, dist_used);

    // If no back references found, fail
    if (dist_used < 2){
        return false;
    }

    // Now get code lengths from the frequencies
    std::vector<u32> ll_nz_counts(ll_used, 0);
    std::vector<u32> ll_nz_symbols(ll_used, 0);
    std::vector<u32> dist_nz_counts(dist_used, 0);
    std::vector<u32> dist_nz_symbols(dist_used, 0);
    u32 ll_used_2 = 0;
    u32 dist_used_2 = 0;
    for (int i = 0; i < 288; i++){
        if (ll_counts.at(i) > 0){
            ll_nz_counts.at(ll_used_2) = ll_counts.at(i);
            ll_nz_symbols.at(ll_used_2++) = i;
        }
    }
    for (int i = 0; i < 32; i++){
        if (dist_counts.at(i) > 0){
            dist_nz_counts.at(dist_used_2) = dist_counts.at(i);
            dist_nz_symbols.at(dist_used_2++) = i;
        }
    }
    assert(ll_used_2 == ll_used);

    if (ll_nz_symbols.at(ll_nz_symbols.size() - 1) < 260 || dist_nz_symbols.at(dist_nz_symbols.size() - 1) < 2)
        return false;


    // Given the (nonzero) frequencies, get the code lengths
    calc_huff_lens(ll_nz_counts, ll_nz_counts.size());
    calc_huff_lens(dist_nz_counts, dist_nz_counts.size());

    // If any code lengths are greater than 15, just use a different block type lol
    for (int i = 0; i < ll_nz_counts.size(); i++){
        if (ll_nz_counts.at(i) > 14){
            return false;
        }
    }

    // Put back into vector of original sizes
    std::vector<u32> ll_code_lengths(288, 0);
    std::vector<u32> dist_code_lengths(32, 0);
    for (int k = 0; k < ll_code_lengths.size(); k++){
        std::vector<u32>::iterator it = find (ll_nz_symbols.begin(), ll_nz_symbols.end(), k);
        if (it != ll_nz_symbols.end()){
            u32 index = distance(ll_nz_symbols.begin(), it);
            u32 code_length = ll_nz_counts.at(index);
            ll_code_lengths.at(k) = code_length;
        }else{
            ll_code_lengths.at(k) = 0;
        }
    }

    for (int k = 0; k < dist_code_lengths.size(); k++){
        std::vector<u32>::iterator it = find (dist_nz_symbols.begin(), dist_nz_symbols.end(), k);
        if (it != dist_nz_symbols.end()){
            u32 index = distance(dist_nz_symbols.begin(), it);
            u32 code_length = dist_nz_counts.at(index);
            dist_code_lengths.at(k) = code_length;
        }else{
            dist_code_lengths.at(k) = 0;
        }
    }

    // No failures --> use type 2!
    stream.push_bit(is_last ? 1 : 0);
    stream.push_bits(2, 2);

    // Write the CL stuff
    size_t size = ll_nz_symbols.at(ll_nz_symbols.size() - 1);
    ll_code_lengths.resize(size + 1);
    size = dist_nz_symbols.at(dist_nz_symbols.size() - 1);
    dist_code_lengths.resize(size + 1);
    write_cl_data(stream, ll_code_lengths, dist_code_lengths);

    std::vector<u32> ll_code = construct_canonical_code(ll_code_lengths);
    std::vector<u32> dist_code = construct_canonical_code(dist_code_lengths);
    type2_LZSS(stream, block_data, block_size, ll_code, ll_code_lengths, dist_code, dist_code_lengths);
    return true;
}

void init(){
    // Initalize the fixed LL code for block type 1
    // Symbols 0 - 143 have length 8
    // Symbols 144 - 255 have length 9
    // Symbols 256 - 279 have length 7
    // Symbols 280 - 287 have length 8
    u32 i;
    for (i = 0; i < 144; i++)
        type1_ll_lengths.at(i) = 8;
    for (i = 144; i < 256; i++)
        type1_ll_lengths.at(i) = 9;
    for (i = 256; i < 280; i++)
        type1_ll_lengths.at(i) = 7;
    for (i = 280; i < 288; i++)
        type1_ll_lengths.at(i) = 8;
    type1_ll_code = construct_canonical_code(type1_ll_lengths);

    // Initialize the fixed distance code for block type 1
    // All distances are represented by 5-bit codes
    for (i = 0; i < 32; i++)
        type1_distance_lengths.at(i) = 5;
    type1_distance_code = construct_canonical_code(type1_distance_lengths);

    // Initialize the base lengths for block type 1
    // If base length <= 10, there are 0 offset bits
    base_length_offset_bits.at(0) = 10;
    // 11 <= length <= 18 there is 1 offset bit --> offset = (length - 11) % 2
    base_length_offset_bits.at(1) = 18;
    // 19 <= length <= 34 there are two offset bits --> offset = (length - 19) % 4
    base_length_offset_bits.at(2) = 34;
    // 35 <= length <= 66 there are three offset bits --> offset = (length - 35) % 8
    base_length_offset_bits.at(3) = 66;
    // 67 <= length <= 130 there are four offset bits --> offset = (length - 67) % 16
    base_length_offset_bits.at(4) = 130;
    // 131 <= length <= 257 there are five offset bits --> offset = (length - 131) % 32
    base_length_offset_bits.at(5) = 257;
    // If length = 258 there are 0 offset bits

    // Initialize the base length symbols so that it's easy to get symbol number from base length
    for (i = 0; i < 9; i++)
        base_length_symbol_nums.at(i) = i + 3;
    for (i = 9; i < 13; i++)
        base_length_symbol_nums.at(i) = i + 3 + (i-8);
    for (i = 13; i < 17; i++)
        base_length_symbol_nums.at(i) = 23 + (i-13)*4;
    for (i = 17; i < 21; i++)
        base_length_symbol_nums.at(i) = 43 + (i-17)*8;
    for (i = 21; i < 25; i++)
        base_length_symbol_nums.at(i) = 83 + (i-21)*16;
    for (i = 25; i < 28; i++)
        base_length_symbol_nums.at(i) = 163 + (i-25)*32;

    // Initialize the base distance symbols
    for (i = 0; i < 5; i++)
        base_dist_symbol_nums.at(i) = i + 1;
    for (i = 6; i < 29; i += 2)
        base_dist_symbol_nums.at(i) = (1 << (i /2)) + 1;
    for (i = 5; i < 30; i += 2)
        base_dist_symbol_nums.at(i) = ((1 << (i / 2)) * 1.5) + 1;
}

int main(){

    init();
    //See output_stream.hpp for a description of the OutputBitStream class
    OutputBitStream stream {std::cout};

    //Pre-cache the CRC table
    auto crc_table = CRC::CRC_32().MakeTable();

    //Push a basic gzip header
    stream.push_bytes( 0x1f, 0x8b, //Magic Number
        0x08, //Compression (0x08 = DEFLATE)
        0x00, //Flags
        0x00, 0x00, 0x00, 0x00, //MTIME (little endian)
        0x00, //Extra flags
        0x03 //OS (Linux)
    );


    //This starter implementation writes a series of blocks with type 0 (store only)
    //Each store-only block can contain up to 2**16 - 1 bytes of data.
    //(This limit does NOT apply to block types 1 and 2)
    //Since we have to keep track of how big each block is (and whether any more blocks 
    //follow it), we have to save up the data for each block in an array before writing it.
    

    //Note that the types u8, u16 and u32 are defined in the output_stream.hpp header
    Block block_contents {};
    u32 block_size {0};
    u32 bytes_read {0};

    char next_byte {}; //Note that we have to use a (signed) char here for compatibility with istream::get()

    //We need to see ahead of the stream by one character (e.g. to know, once we fill up a block,
    //whether there are more blocks coming), so at each step, next_byte will be the next byte from the stream
    //that is NOT in a block.

    //Keep a running CRC of the data we read.
    u32 crc {};


    if (!std::cin.get(next_byte)){
        //Empty input?
        
    }else{

        bytes_read++;
        //Update the CRC as we read each byte (there are faster ways to do this)
        crc = CRC::Calculate(&next_byte, 1, crc_table); //This call creates the initial CRC value from the first byte read.
        //Read through the input
        while(1){
            block_contents.at(block_size++) = next_byte;
            if (!std::cin.get(next_byte))
                break;

            bytes_read++;
            crc = CRC::Calculate(&next_byte,1, crc_table, crc); //Add the character we just read to the CRC (even though it is not in a block yet)

            //If we get to this point, we just added a byte to the block AND there is at least one more byte in the input waiting to be written.
            if (block_size == block_contents.size()){
                //The block is full, so write it out.
                //We know that there are more bytes left, so this is not the last block
                bool ret = generate_type2(stream, block_contents, block_size, false);
                if (! ret){
                    generate_type1(stream, block_contents, block_size, false);
                }
                block_size = 0;
            }
        }
    }

    //At this point, we've finished reading the input (no new characters remain), and we may have an incomplete block to write.
    if (block_size > 0){
        //Write out any leftover data
        bool ret = generate_type2(stream, block_contents, block_size, true);
        if (!ret){
            generate_type1(stream, block_contents, block_size, true);
        }
        block_size = 0;
    }
    stream.flush_to_byte();

    //Now close out the bitstream by writing the CRC and the total number of bytes stored.
    stream.push_u32(crc);
    stream.push_u32(bytes_read);

    return 0;
}