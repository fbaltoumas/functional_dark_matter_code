from sys import stdout, stderr
import multiprocessing as mp
import argparse as ap
import subprocess as sp
import time

# TMalign and MMalign locations
tmalign = "/usr/local/bin/TMalign"
mmalign = "/opt/mmalign/MMalign"

# GLOBALS --- these can be overriden by the parser arguments
method = tmalign
cpus = mp.cpu_count()


def create_parser():
    parser = ap.ArgumentParser(description="Compare two lists of PDB files with TMalign/MMalign in parallel.",
                               add_help=False,
                               formatter_class=ap.ArgumentDefaultsHelpFormatter)
    required = parser.add_argument_group('required arguments')
    required.add_argument("--query", type=str, action="store", help="list of query PDB files", required=True)
    required.add_argument("--target", type=str, action="store", help="list of target PDB files", required=True)

    optional = parser.add_argument_group('optional arguments')
    optional.add_argument('--method', type=str, action='store', help='alignment method: tmalign or mmalign',
                          default='tmalign', required=False)
    optional.add_argument('--cpus', type=int, action='store', help='No. of CPUs to use. Set to 0 to use all CPUs',
                          default=0, required=False)
    optional.add_argument('-h', '--help', action='help', help='Show this help message.')
    args = parser.parse_args()
    return args


def parse_list(inpfile):
    fl = open(inpfile, "r")
    fl_lines = [i.rstrip() for i in fl]
    fl.close()
    return fl_lines


def align_to_list(query_pdb, target_list, method):
    # TMalign and MMalign outputs differ by ONE word in the lines containing PDB lengths:
    # in TMalign, it says "Chain", in MMalign, it says "Structure"
    if method == tmalign:
        method_prefix = "Chain"
    else:
        method_prefix = "Structure"
    for target_pdb in target_list:
        # format the shell command to call
        cmd = "%s %s %s -a T" % (method, query_pdb, target_pdb)
        # call the shell command and assign it to proc
        proc = sp.Popen(cmd, shell=True, stdout=sp.PIPE, stderr=stderr)
        # parse the output
        out = proc.communicate()
        # get the STDOUT of the command and decode it to UTF-8. Also, split it based on newlines
        result = out[0].decode('utf-8').split("\n")
        # print the output (containing the alignment itself) in the STDERR
        print("\n".join(result), file=stderr)
        # initializing the TM-scores list
        scores = []
        # parsing the output
        for r in result:
            # detect lines containing TM-score values and append them
            if r != "" and r.split()[0] == "TM-score=":
                scores.append(float(r.split()[1]))
            # detect the line containing the length of query_pdb
            if r != "" and r.split(":")[0] == "Length of %s_1" % method_prefix:
                c1 = r.split(":")[1].split()[0]
            # detect the line containing the length of target_pdb
            if r != "" and r.split(":")[0] == "Length of %s_2" % method_prefix:
                c2 = r.split(":")[1].split()[0]
        if len(scores) > 0:
            res = [query_pdb, target_pdb, str(c1), str(c2), "%.3f" % scores[0], "%.3f" % scores[1], "%.3f" % scores[2]]
            # use stdout.write instead of print to avoid missing newlines during threading (race conditions etc.)
            stdout.write("\t".join(res) + " ;\n")


if __name__ == "__main__":
    start_time = time.perf_counter()
    arg_parser = create_parser()  # create the arguments parser
    # parse CMD input
    query_list = parse_list(arg_parser.query)
    target_list = parse_list(arg_parser.target)
    if arg_parser.method.lower() == "tmalign":
        method = tmalign
    else:
        method = mmalign
    if arg_parser.cpus > 0:
        cpus = arg_parser.cpus
    # pool resources for parallel execution, based on the assigned cpu count
    pool = mp.Pool(cpus)
    # distribute jobs among pooled resources asynchronously
    proc = [pool.apply_async(align_to_list, args=(x, target_list, method,)) for x in query_list]
    # get the job results
    result = [p.get() for p in proc]
    end_time = time.perf_counter()
    print("Total Finished in %f seconds" % (end_time - start_time + 1), file=stderr)


