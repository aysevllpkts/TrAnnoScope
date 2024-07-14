import sys

def reformat_header(header):
    parts = header.split('; ')
    transcript_id = parts[0].split()[0][1:]
    length_info = parts[1].split('=')[1].split(',')[0]
    orf_info = parts[1].split(',')[2]
    coords = parts[3].split('=')[1]
    strand = "+"                        # evigene does not have strand info for all reads, I put (+) for all of them.
    score = parts[2].split('=')[1]      # just to mimick transdecoder header, it is not the same thing. My version shows the clen as score.

    new_header = f">{transcript_id}.p1 GENE.{transcript_id}~~{transcript_id}.p1 ORF type:{orf_info} len:{length_info} ({strand}),score={score} {transcript_id}:{coords}({strand})\n"

    return new_header

def reformat_fasta(input_file, output_file):
    with open(input_file, 'r') as f_in, open(output_file, 'w') as f_out:
        for line in f_in:
            if line.startswith('>'):
                new_header = reformat_header(line.strip())
                f_out.write(new_header)
            else:
                f_out.write(line)


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python script.py input_file output_file")
        sys.exit(1)

    input_file = sys.argv[1]
    output_file = sys.argv[2]
    reformat_fasta(input_file, output_file)
    print("Headers reformatted and saved to output file.")


#reformat_fasta(snakemake.input[0], snakemake.output[0])
#print("Headers reformatted and saved to output file.")
