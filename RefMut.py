import sys, re, subprocess, copy, argparse
from optparse import OptionParser



class RefMut(object):

    def __init__(self, fasta, max_mut=5, incr_mut=1):

        self.fasta = fasta
        self.directory_path = ""
        self.file_root = ""
        self.ref_seq = ""
        self.header = ""
        self.ref_seq_len = 0

        self.read_fasta()
        self.create_environment()

        if 0 < max_mut < 100:
            self.max_mut = max_mut
        else:
            print "The maximum rate of mutation must be between 0 and 100.\n Exiting RefMut..."
            sys.exit()

        if incr_mut >= 0 and (max_mut % incr_mut) == 0:
            self.incr_mut = incr_mut
        else:
            print "The increment of mutation rate must be 0 or greater and be a factor of the maximum rate of mutation." \
                  "\nExiting RefMut..."
            sys.exit()

        self.mutate()

    def read_fasta(self):
        """ Read Fasta file and log the length of the reference sequence"""

        ffh = open(self.fasta, 'r')
        root = re.match('(.+)(.fasta)', self.fasta)
        self.file_root = root.group(1)
        self.header = ffh.readline().strip()
        for line in ffh:
            self.ref_seq += line.strip()
        self.ref_seq_len = len(self.ref_seq)
        ffh.close()

        return
    
    def create_environment(self):
        """ Create a directory for mutated file outputs """
        
        self.directory_path = "RefMut_" + self.file_root
        print self.directory_path
        cmd = ['mkdir', '-p', self.directory_path]
        subprocess.call(cmd)

        return
        
    def mutate(self):
        """ loop through mutation rates by increment rate until max mutation rate is reached """
        
        mutations = ['A', 'T', 'G', 'C']
        mutated_positions = []
        mutations_number_log = [0]

        for mut_rate in range(0, self.max_mut + 1, self.incr_mut):
            # create new file for this mutated reference sequence
            filename = '%s/%s_%sperc.fasta' % (self.directory_path, self.file_root, mut_rate)
            ofh = open(filename, 'w')
            number_of_mutations = mut_rate * self.ref_seq_len / 100
            new_mutations = number_of_mutations - mutations_number_log[len(mutations_number_log) - 1]
            out = 'Mutated %i%% with %i total mutations and %i new mutations' % (mut_rate, number_of_mutations, new_mutations)
            print out
            for mutation_instance in range(new_mutations):
                position = 0
                new_pos = False
                while new_pos is False:
                    position = randint(0, self.ref_seq_len - 1)
                    if position not in mutated_positions:
                        mutated_positions.append(position)
                        new_pos = True
                #print self.ref_seq_len, position
                current_symbol = self.ref_seq[position]
                temp_table = copy.copy(mutations)
                temp_table.remove(current_symbol)
                mutation = temp_table[randint(0, 2)]
                self.ref_seq = self.ref_seq[0:position] + mutation + self.ref_seq[position +1:]
                #print "Position", position, current_symbol, ">", mutation

            new_header = '%s, %i PERCENT MUTATED' % (self.header, mut_rate)
            print >> ofh, new_header
            print >> ofh, self.ref_seq
            ofh.close()
            # add number of mutations to log to be reference during next loop
            mutations_number_log.append(number_of_mutations)
        return


def main():
    """ 
    Command line method if script is implemented from the commandline and not imported as a class.
    """
    
    parser = argparse.ArgumentParser(description='RefMut mutates a given reference into new fasta files that are incrementally different from the originial given sequence. ')
    parser.add_argument('-R', dest='refseq', help='Fasta reference sequence to be mutated', required=True)
    parser.add_argument('-M', dest='maxi', help='Maximum rate of mutation', type=int, default=5)
    parser.add_argument('-I', dest='incr', help='Increment of mutation rate', type=int, default=1)
    results = parser.parse_args()
    
    # assign command line args to variables
    reference_sequence = results.refseq
    max_mutation = results.maxi
    incr_mutation = results.incr
    print reference_sequence, max_mutation, incr_mutation
    RefMut(reference_sequence, max_mutation, incr_mutation)

if __name__ == '__main__':
    main()
#    RefMut('EcoliK12_refseq.fasta', 20, 5)
