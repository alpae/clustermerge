
import gzip
import os
import argparse
import re
import collections
import lark
import csv
import itertools
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt



class ListComparer:
    # Comparing two lists of elements by casting them to a set => each element occurs only once
    def __init__(self, list1, list2):
        self.list1 = list1
        self.list2 = list2
        self.set1 = set(list1)
        self.set2 = set(list2)

    def get_union(self):
        return list(self.set1 | self.set2)

    def get_intersection(self):
        return list(self.set1 & self.set2)

    def only_in_list1(self):
        return list(self.set1 - self.set2)

    def only_in_list2(self):
        return list(self.set2 - self.set1)

    def get_total_number(self):
        return len(self.get_union())

    def get_total_number_in_both(self):
        return len(self.get_intersection())

    def get_total_number_in_list1(self):
        return len(self.set1)

    def get_total_number_in_list2(self):
        return len(self.set2)

    def list1_has_duplicates(self):
        return len(self.get_duplicates_of_list1()) > 0

    def list2_has_duplicates(self):
        return len(self.get_duplicates_of_list2()) > 0

    def get_duplicates_of_list1(self):
        return list(set([x for x in self.list1 if self.list1.count(x) > 1]))

    def get_duplicates_of_list2(self):
        return list(set([x for x in self.list2 if self.list2.count(x) > 1]))

Pair = collections.namedtuple('Pair', ['Entry1', 'Entry2', 'Rng1', 'Rng2', 'Score'])

grammar = '''ref: ( ref_fun )*
             ?ref_fun: "RefinedMatches(" matches ")" _TERM
             matches: "[" [ match ("," match)* ] "]"
             match : "NULL"
                   | "[" int_ "," int_ "," real_ "," real_ "," rng "," rng "," real_ [ "," int_ ] "]"
             ?int_  : /[0-9]+/     -> int_
             ?real_ : /[0-9.]+/    -> real
             rng    : int_ _RNG int_
             _RNG   : ".."
             _TERM : /[:;]/
             %import common.WS
             %ignore WS
             '''


class MatchTrans(lark.Transformer):
    def int_(self, vals):
        return int(vals[0])

    def rng(self, vals):
        return vals[0], vals[1]

    def real(self, vals):
        return float(vals[0])

    def match(self, vals):
        if len(vals) >= 5:
            return Pair(vals[0], vals[1], vals[4], vals[5], int(vals[2]))

    def matches(self, vals):
        return [v for v in vals if v is not None]

    def ref(self, matches):
        return itertools.chain.from_iterable(matches)


def process_file(filename):
    print("file is {}".format(filename))
    open_ = gzip.open if filename.endswith('.gz') else open
    parser = lark.Lark(grammar, start='ref', parser='lalr', transformer=MatchTrans())
    with open_(filename) as fh:
        data = "\n".join([line for line in fh if not line.startswith('#')])
        matches = parser.parse(data)
    return matches


def load_matches(folderpath, minscore=181, genome_pairs=None):
    result = collections.defaultdict(list)
    folders = [x[0] for x in os.walk(folderpath) if x[0] != folderpath]
    for folder in folders:
        genome1 = folder.split("/")[-1]
        files = [f for f in os.listdir(folder) 
                 if os.path.isfile(os.path.join(folder, f)) and not '.sha2.' in f]
        for fname in files:
            genome2 = fname.split("/")[-1].split("_")[0].split(".gz")[0]
            full_file_path = os.path.join(folder, fname)
            if genome_pairs is not None and (genome1, genome2) not in genome_pairs:
                 print("skipping genome pair file {}".format(full_file_path))
                 continue
            candidates = [pair for pair in process_file(full_file_path) if pair.Score >= minscore]
            key = (genome1, genome2)
            result[key].extend(candidates)
    return result


def main(protclusterdir, allalldir, minscore):
    data_dict = load_matches(protclusterdir, minscore)
    ref_dict = load_matches(allalldir, minscore, data_dict.keys()) 

    if set(ref_dict) != set(data_dict):
        print(set(ref_dict))
        print(set(data_dict))
        print("Lengths: {} and {}".format(len(set(ref_dict)), len(set(data_dict))))
        print(set(ref_dict).symmetric_difference(set(data_dict)))
        raise ValueError("The dictionaries do not have the same keys (genomes)!")

    reported_by_both = 0
    reported_by_ref_only = 0
    reported_by_data_only = 0
    total_pairs_in_ref = 0
    total_pairs_in_data = 0

    additional_pairs_in_ref = dict()
    additional_pairs_in_data = dict()

    for ref_key in ref_dict:
        list_comp = ListComparer(ref_dict[ref_key], data_dict[ref_key])

        reported_by_both += list_comp.get_total_number_in_both()
        reported_by_ref_only += len(list_comp.only_in_list1())
        reported_by_data_only += len(list_comp.only_in_list2())
        total_pairs_in_ref += list_comp.get_total_number_in_list1()
        total_pairs_in_data += list_comp.get_total_number_in_list2()

        additional_pairs_in_ref[ref_key] = list_comp.only_in_list1()
        additional_pairs_in_data[ref_key] = list_comp.only_in_list2()
        for dubious in list_comp.only_in_list2():
            for possible_cand in list_comp.only_in_list1():
                if dubious.Entry1 == possible_cand.Entry1 and dubious.Entry2 == possible_cand.Entry2:
                    print("match only in data found with modification in ref:\n  {} vs {}".format(dubious, possible_cand))
                    break
                

    print("additional in data:\n {}".format(additional_pairs_in_data))
    print("additional in ref:\n {}".format(additional_pairs_in_ref))
    
    
    missed_scores = []
    dataset_name = protclusterdir.split('/')[-1]
    with open(os.path.join(protclusterdir, 'missed.txt'), 'w') as fout:
        csv_writer = csv.writer(fout, delimiter='\t')
        csv_writer.writerow(["Dataset","Genome1","Genome2","Entry1","Entry2","Score"])

        for key, value in additional_pairs_in_ref.iteritems():
            #print("key {}, value {}".format(key, value))
            # iterate over values, get all missed scores
            for v in value:
                missed_scores.append(int(v[-1]))
                csv_writer.writerow([dataset_name, key[0], key[1], v.Entry1, v.Entry2, v.Score])

    missed_scores.sort(reverse=True)
    print("missed scores are {}".format(missed_scores))
    fig = plt.hist(missed_scores, 'scott')
    plt.title("ClusterMerge: Scores of Missed Matches")
    plt.xlabel("Score")
    plt.ylabel("Number of Missed Matches")
    plt.savefig("hist-{}.png".format(protclusterdir.split('/')[-1]))
    print("average missed score is {}".format(sum(missed_scores) / len(missed_scores)))
    print("median missed score is {}".format(missed_scores[len(missed_scores)/2]))

    percent_both = round(float(reported_by_both) / total_pairs_in_ref, 7) * 100
    percent_reference = round(float(reported_by_ref_only) / total_pairs_in_ref, 7) * 100
    percent_data = round(float(reported_by_data_only) / total_pairs_in_ref, 7) * 100

    output = ""
    output += "Reported by both: {} ({}%)\n".format(reported_by_both, percent_both)
    output += "Reported by reference only: {} ({}%)\n".format(
        reported_by_ref_only, percent_reference
    )
    output += "Reported by new only: {} ({}%)\n".format(
        reported_by_data_only, percent_data
    )
    output += "Total reported significands: {}".format(total_pairs_in_data)
    print("total in ref {}".format(total_pairs_in_ref))
    print("total in data {}".format(total_pairs_in_data))

    print(output)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Compare two Darwin result DBs.")
    parser.add_argument('--minscore','-s', default=181, type=float, help="Min score threshold to filter allall files")
    parser.add_argument("allalldir", help="The dir with gz files from allall")
    parser.add_argument("protclusterdir", help="The dir with files from protcluster")

    conf = parser.parse_args()

    main(conf.protclusterdir, conf.allalldir, conf.minscore)
