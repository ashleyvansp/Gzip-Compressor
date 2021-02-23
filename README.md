# DEFLATE
ASSIGNMENT 2 - CSC485B
Ashley Van Spankeren


A summary of the functions you'll find in the fantastic 900+ lines of code:
* construct_canonical_code: I borrowed this directly from the start package provided. Unmodified.
* write_ll_code_and_offset: This is used when a backreference is found; it outputs the base length and offset directly to the stream.
* write_dist_code_and_offset: Same as above; outputs the base distance and offset (if applicable) directly to the stream.
* generate_type1: Handles everything to convert block_contents into a valid type 1 block, including running LZSS.
* type2_LZSS_frequencies: Called in the initial phases of writing a type 2 block. Runs LZSS on the block contents and updates a vector of the frequencies of the lengths/ literals/ distances found.
* get_dist_freq: Gets called from type2_LZSS_frequencies; finds the base distance of the 'dist' passed in and increments its frequency in dist_counts.
* get ll_freq: Similar to above; finds the base length of the 'len' passed in and increments its frequency in ll_counts.
* type2_LZSS: Runs LZSS on the block_data and outputs directly to stream using codes from the ll_code and dist_code.
* KM_inequality: Returns true if and only if the Kraft-McMillan value of the lengths in 'W' equals to 1.
* make_optimal: Given a vector of code lengths, this function essentially 'rearranges the tree' of lengths until they satisfy the Kraft-McMillan inequality.
* calc_huff_lens: Takes a vector of frequencies as input and operates in place to convert each frequency to its optimal code length so that the Kraft-McMillan inequality is satisfied.
* rle: Applies run-length encoding to a vector of code lengths so that they can be encoded with the CL code (ie. collapes runs of zeroes into the symbol 18 & the run length).
* write_code_length_encoding: Takes the RLE vector from above, the CL codes, and the CL code lengths, and writes to stream the encoding of the code lengths.
* get_cl_freqs: Given the RLE of the LL and distance code lengths, this function updates the vector based on the frequency of each CL symbol.
* write_cl_data: Modified from the starter code provided. Calculates HLIT, HDIST, and HCLEN, and makes a CL code accordingly, then outputs it all to stream.
* generate_type2: Handles everything to convert block_data into a valid type 2 block.
* init: Initializes the bases values of the type 1 LL and distance code.
* main: Reads from the stream of input and decides the block type to use.

Features implemented:
* This compressor uses both block types 1 and 2 to compress the data. 
* It first tries to write the data using block type 2, but contains logic so that if 'red flags' are encountered (such as very few backreferences, or very short backreferences) then it reverts to type 1.
* Uses a window size of 2500 for searching for backreferences. This can easily be changed by modifying the global variable window_size.
* Uses a fixed block length of BLOCK_MAX.
* The time to compress the entire collection of test data is around 15 seconds (when ran on the school server).
* Generates dynamic LL and distance codes (in the function generate_type2) by:
    - Running LZSS on the block data
    - Keeping track of the frequencies of each literal, length, and distance symbol
    - Generating a customized Huffman code so that only the LL and distance symbols used in this block are allocated a code
* Optimizes the block type 2 header (in the function write_cl_data) by:
    - Running RLE on the LL and distance code length vectors
    - Analyzing which CL codes are used in the run length encodings
    - Generating a customized CL Huffman code so that only the CL symbols used in this block's RLE are allocated a code
* Optimizes the choice of block type (in the function generate_type2, and in main) by avoiding block type 2 (and its lengthy header) if the LZSS frequency vector shows that few backreferences are used (or that the backreferences have short lengths)

Admittedly I didn't use fancy data structures in this project, so there's little to write about in that department.

Notes about some algorithms used:
* I chose a very naive approach to LZSS:
    - beginning at the start of the sliding window, look for backreferences
    - as soon as one is found, obtain the base length, length offset, base distance, and distance offset and write them to the stream
    - move to the next character, and repeat
* The algorithm I used to generate the lengths of the Huffman code was taken directly from the pseudocode from Alistair Moffat's 'Huffman Coding' (see reference below)
* The run length encoding algorithm is quite simple; it collects the longest possible run starting from the current index then collapses the run into a CL code (if run longer than 2) or outputs the literal otherwise

Things I would've liked to include, but ran out of time for:
* As mentioned, LZSS finds the first available backreference in the sliding window and outputs it immediately; I would've liked for the algorithm to continue looking for backreferences with longer length or shorter distance
* Rather than outputting the block directly to the stream when running LZSS, I would've preferred to store it somehow so that the bitstreams for type 1 and type 2 could be compared and the shorter one gets chosen to be written to the stream.
* Minor bug fixes, such as the CL code lengths sometimes surpass 7 and instead of dealing with this like a real adult I chose instead to revert to a fixed code length.
* Adding in logic so that block type 0 is chosen if it's most efficient for a particular block
* Implenting the Package-Merge algorithm to generate the Huffman code
* General refactoring....

Sources used:
* https://thispointer.com/c-how-to-find-an-element-in-vector-and-get-its-index/ - to search for a value in a vector and obtain its index
* https://people.eng.unimelb.edu.au/ammoffat/abstracts/compsurv19moffat.pdf - to calculate the Huffman code lengths, given the frequencies
* https://tools.ietf.org/html/rfc1951#section-3.2.7 - for general reference and clarification
