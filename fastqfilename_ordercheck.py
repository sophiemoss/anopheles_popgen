import subprocess

with open('samples.txt', 'r') as file:
    lines = file.readlines()

    for line in lines:
        SID = line.strip()
        print (SID)
        output_file_1 = f"{SID}_read1_names.txt"
        output_file_2 = f"{SID}_read2_names.txt"

        command1 = f"seqkit seq --name {SID}_1.fastq.gz | awk '{{print $1}}' > {output_file_1}"
        command2 = f"seqkit seq --name {SID}_2.fastq.gz | awk '{{print $1}}' > {output_file_2}"
        command3 = f"diff {output_file_1} {output_file_2}  | wc -l > {SID}_difference.log"

        subprocess.run(command1, shell=True)
        subprocess.run(command2, shell=True)
        subprocess.run(command3, shell=True)

        countlines = open(f'{SID}_difference.log', 'r').readline().strip()

        print (SID,countlines)


