import msprime
import tskit
import tqdm
import numpy as np
from IPython.core.display import HTML
from jupyterquiz import display_quiz


class DownloadProgressBar(tqdm.tqdm):
    def update_to(self, b=1, bsize=1, tsize=None):
        if tsize is not None:
            self.total = tsize
        self.update(b * bsize - self.n)


class Workbook:
    css = """<style>
        dl {border: green 1px solid; margin-top: 1em}
        dt {color: white; background-color: green; padding: 4px; display: block; }
        dd {padding: 4px;}
    </style>"""

    # Used for making SVG formatting smaller
    small_class = "x-lab-sml"
    small_style = (
        ".x-lab-sml .sym {transform:scale(0.6)} .x-lab-sml .lab {font-size:7pt;}" # All labels small
        ".x-lab-sml .x-axis .tick .lab {"
        "font-weight:normal;transform:rotate(90deg);text-anchor:start;dominant-baseline:central;}"
    )

    # some useful functions
    @staticmethod
    def convert_metadata_to_new_format(ts):
        # Quick hack to read individual and population metadata as a python dict
        tables = ts.dump_tables()
        tables.populations.metadata_schema = tskit.MetadataSchema.permissive_json()
        tables.individuals.metadata_schema = tskit.MetadataSchema.permissive_json()
        tables.sites.metadata_schema = tskit.MetadataSchema.permissive_json()
        return tables.tree_sequence()

    @staticmethod
    def download(url):
        return DownloadProgressBar(
            unit='B', unit_scale=True, miniters=1, desc=url.split('/')[-1])



class Workbook1(Workbook):

    def __init__(self):
        self.ts = self.simulate_ts()
        assert len(self.ts.site(12).mutations) == 2  # check there are sites with multiple mutations
        self.ts.dump("simulated.trees")

    def simulate_ts(self):
        pop_size=1000
        seq_length=1_000_000
        
        sweep_model = msprime.SweepGenicSelection(
            position=seq_length/2, start_frequency=0.001, end_frequency=0.999, s=0.25, dt=1e-6)
        
        ts = msprime.sim_ancestry(
            10,
            model=[sweep_model, msprime.StandardCoalescent()],
            population_size=pop_size,
            sequence_length=seq_length,
            recombination_rate=1e-8,
            random_seed=654321,  # only needed for repeatabilty
        )
        ts = msprime.sim_mutations(ts, rate=2e-8, random_seed=203)
        tables = ts.dump_tables()
        tables.individuals.metadata_schema = tskit.MetadataSchema.permissive_json()
        tables.individuals.packset_metadata([
            tables.individuals.metadata_schema.validate_and_encode_row({"name": n})
            for n in ["Ada", "Bob", "Cat", "Dee", "Eli", "Fi", "Guy", "Hal", "Ida", "Jo"]
        ])
        return tables.tree_sequence()
        
    def Q1(self):
        display_quiz([
            {
                "question": "How many edges in this tree sequence?",
                "type": "numeric",
                "answers": [
                    {
                        "type": "value",
                        "value": self.ts.num_edges,
                        "correct": True,
                        "feedback": "Correct."
                    },
                    {
                        "type": "range",
                        "range": [ -100000000, 100000], 
                        "correct": False,
                        "feedback": "Try again (hint: look at `ts.num_edges()`)"
                    },
                ]
            },
            {
                "question": "How many sites in this tree sequence?",
                "type": "numeric",
                "answers": [
                    {
                        "type": "value",
                        "value": self.ts.num_sites,
                        "correct": True,
                        "feedback": "Correct."
                    },
                    {
                        "type": "range",
                        "range": [ -100000000, 100000], 
                        "correct": False,
                        "feedback": "Try again (hint: look at `ts.num_sites()`)"
                    },
                ]
            },
            {
                "question": "How many mutations in this tree sequence?",
                "type": "numeric",
                "answers": [
                    {
                        "type": "value",
                        "value": self.ts.num_mutations,
                        "correct": True,
                        "feedback": 
                            "Correct: there may be more mutations than sites because "
                            "there could be multiple mutations at a single site"
                    },
                    {
                        "type": "range",
                        "range": [ -100000000, 100000], 
                        "correct": False,
                        "feedback": "Try again (hint: look at `ts.num_mutations()`)"
                    },
                ]
            },
        ])

    def Q2(self):
        display_quiz([
            {
                "question": "How many mutations at site 11?",
                "type": "numeric",
                "answers": [
                    {
                        "type": "value",
                        "value": len(self.ts.site(11).mutations),
                        "correct": True,
                        "feedback": "Correct."
                    },
                    {
                        "type": "range",
                        "range": [ -100000000, 100000], 
                        "correct": False,
                        "feedback": "Try again (hint: look at `ts.site(11)`)"
                    },
                ]
            },
            {
                "question": "How many mutations at site 12:",
                "type": "numeric",
                "answers": [
                    {
                        "type": "value",
                        "value": len(self.ts.site(12).mutations),
                        "correct": True,
                        "feedback": "Correct."
                    },
                    {
                        "type": "range",
                        "range": [ -100000000, 100000], 
                        "correct": False,
                        "feedback": "Try again (hint: look at `ts.site(12)`)"
                    },
                ]
            },
        ])

    def Q3(self):
        tree = self.ts.at(400_000)
        mut_counts = np.bincount([tree.num_samples(m.node) for m in tree.mutations()])
        display_quiz([
            {
                "question":
                    "How many mutations in total are there in the tree at position "
                    "400 Kb?",
                "type": "numeric",
                "answers": [
                    {
                        "type": "value",
                        "value": tree.num_mutations,
                        "correct": True,
                        "feedback": "Correct."
                    },
                    {
                        "type": "range",
                        "range": [ -1000000, 1000000], 
                        "correct": False,
                        "feedback":
                            "Try again"
                    },
                ]
            },
            {
                "question":
                    "How many singleton mutations (above one sample, i.e. on a terminal "
                    "branch) are there in the tree at position 400 Kb?",
                "type": "numeric",
                "answers": [
                    {
                        "type": "value",
                        "value": int(mut_counts[1]),
                        "correct": True,
                        "feedback": "Correct."
                    },
                    {
                        "type": "range",
                        "range": [ -1000000, 1000000], 
                        "correct": False,
                        "feedback":
                            "Try again"
                    },
                ]
            },
            {
                "question":
                    "How many doubleton mutations (above 2 samples) are there in the "
                    "tree at position 400 Kb?",
                "type": "numeric",
                "answers": [
                    {
                        "type": "value",
                        "value": int(mut_counts[2]),
                        "correct": True,
                        "feedback": "Correct."
                    },
                    {
                        "type": "range",
                        "range": [ -1000000, 1000000], 
                        "correct": False,
                        "feedback":
                            "Try again"
                    },
                ]
            },
        ])

    def Q4(self):
        display_quiz([{
            "question":
                "What is the age (to the nearest generation) of the root in the first "
                "tree?",
            "type": "numeric",
            "precision": 0,
            "answers": [
                {
                    "type": "value",
                    "value": round(self.ts.node(self.ts.first().root).time, 1),
                    "correct": True,
                    "feedback": "Correct."
                },
                {
                    "type": "range",
                    "range": [ -100000000, 1000000], 
                    "correct": False,
                    "feedback":
                        "Try again (hint: the root node in the tree "
                        f"has ID {self.ts.first().root};"
                        " feed that to the `ts.node()` method)"
                },
            ]
        }])

    def Q5a(self):
        correct_name = self.ts.individual(self.ts.node(14).individual).metadata["name"]
        display_quiz([{
            "question":
                "What is the name of the individual associated with sample node 14 "
                "in the original tree sequence?",
            "type": "multiple_choice",
            "answers": [
                {
                    "answer": i.metadata["name"],
                    "correct": i.metadata["name"] == correct_name,
                    "feedback": (
                        "Correct" if i.metadata["name"] == correct_name
                        else "Sorry, that's not right."
                        )
                }
                for i in self.ts.individuals()
            ]
        }])

    def Q5b(self):
        display_quiz([{
            "question":
                "How many populations are defined in this tree sequence?",
            "type": "numeric",
            "precision": 0,
            "answers": [
                {
                    "type": "value",
                    "value": self.ts.num_populations,
                    "correct": True,
                    "feedback":
                        "Correct. All nodes in the tree sequence belong to the same "
                        "population (imaginatively named 'pop_0')"
                },
                {
                    "type": "range",
                    "range": [ -100000000, 1000000], 
                    "correct": False,
                    "feedback":
                        "Try again"
                },
            ]
        }])

    def Q6a(self):
        display_quiz([
            {
                "question":
                    "What alleles does Dee have at the last variable site in the "
                    "simplified tree sequence",
                "type": "multiple_choice",
                "answers": [
                    {"answer": "Homozygous AA", "correct": False, "feedback": "Try again"},
                    {"answer": "Heterozygous AC", "correct": True, "feedback": "Correct"},
                    {"answer": "Homozygous CC", "correct": False, "feedback": "Try again"},
                ]
            },
        ])

    def Q6b(self):
        for i, v in enumerate(self.ts.variants()):
            if i == 0:
                a1 = int(v.genotypes[0])
            if i == self.ts.num_sites-1:
                a2 = int(v.genotypes[1])

        display_quiz([
            {
                "question":
                    "What is the genotypic state of the first sample at the first site?",
                "type": "multiple_choice",
                "answers": [
                    {
                        "answer": f"{s}", "correct": (a1 == s),
                        "feedback": ("Correct" if a1 == s else "Sorry, that's not right.")
                    }
                    for s in [0, 1]
                ]
            },
            {
                "question":
                    "What is the genotypic state of the second sample at the last site?",
                "type": "multiple_choice",
                "answers": [
                    {
                        "answer": f"{s}", "correct": (a2 == s),
                        "feedback": ("Correct" if a2 == s else "Sorry, that's not right.")
                    }
                    for s in [0, 1]
                ]
            },
        ])

class Workbook2(Workbook):
    
    
    def __init__(self):
        small_ts = msprime.sim_ancestry(
            samples=2,
            sequence_length=1e3, # 1Kb
            recombination_rate=1e-8, # allow for recombination
            population_size=20_000, # a rough "effective population size" suitable for humans
            random_seed=107
        )
        self.mts_small = msprime.sim_mutations(small_ts, rate=2e-7, random_seed=103)
    
        # First model
        self.ts1 = msprime.sim_ancestry(
            20_000,
            sequence_length=1e6,
            recombination_rate=1e-8,
            population_size=20_000,
            random_seed=2022,
        )
        self.mts1 = msprime.sim_mutations(self.ts1, rate=1e-8, random_seed=2022)
    
        # Second model
        deme_size = 500 # population size of each deme
        num_demes = 8
        num_deme_samples = 25
        mu = 1e-8
        demography = msprime.Demography.stepping_stone_model(
            [deme_size] * num_demes,
            migration_rate=0.05
        )
        ts = msprime.sim_ancestry(
            {i: num_deme_samples for i in range(num_demes)},
            sequence_length=5e6, # 1 Mbp
            demography=demography,
            recombination_rate=1e-8, # human-like recombination rate
            random_seed=123,
        )
        self.mts =  msprime.sim_mutations(
            ts,
            rate=mu, # human-like mutation rate
            random_seed=321
        )
        
        Fst_values = []
        for ts in msprime.sim_ancestry(
            {i: num_deme_samples for i in range(num_demes)},
            sequence_length=1e6,
            demography=demography,
            recombination_rate=1e-8,
            random_seed=1234,
            num_replicates=100
        ):
            Fst = ts.Fst([ts.samples(0), ts.samples(1)], mode="branch")
            Fst_values.append(float(Fst))

        self.Fst_0_1_mean_100_reps = np.mean(Fst_values)
            
    def Q1(self):
        display_quiz([{
            "question":
                "How many trees are in your newly simulated tree sequence?",
            "type": "numeric",
            "precision": 0,
            "answers": [
                {
                    "type": "value",
                    "value": self.ts1.num_trees,
                    "correct": True,
                    "feedback":
                        "Correct"
                },
                {
                    "type": "range",
                    "range": [ -100000000, 1000000], 
                    "correct": False,
                    "feedback":
                        "Try again"
                },
            ]
        }])
        
    def Q2(self):
        display_quiz([{
            "question":
                "What is the ID of the site with two mutations?",
            "type": "numeric",
            "precision": 0,
            "answers": [
                {
                    "type": "value",
                    "value": [
                        s.id for s in self.mts_small.sites() if len(s.mutations)==2
                    ][0],
                    "correct": True,
                    "feedback":
                        "Correct"
                },
                {
                    "type": "range",
                    "range": [ -100000000, 1000000], 
                    "correct": False,
                    "feedback":
                        "Try again"
                },
            ]
        }])

    def Q3(self):
        display_quiz([
            {
                "question":
                    "How many variable sites are in the tree sequence?",
                "type": "numeric",
                "precision": 0,
                "answers": [
                    {
                        "type": "value",
                        "value": self.mts1.num_sites,
                        "correct": True,
                        "feedback":
                            "Correct"
                    },
                    {
                        "type": "range",
                        "range": [ -100000000, 1000000], 
                        "correct": False,
                        "feedback":
                            "Try again"
                    },
                ]
            },
            {
                "question":
                    "How many big is the tree sequence, in MiB (binary megabytes) "
                    "to 1 decimal place?",
                "type": "numeric",
                "precision": 1,
                "answers": [
                    {
                        "type": "value",
                        "value": round(self.mts1.nbytes/1024/1024, 1),
                        "correct": True,
                        "feedback":
                            "Correct"
                    },
                    {
                        "type": "range",
                        "range": [ -100000000, 1000000], 
                        "correct": False,
                        "feedback":
                            "Try again (hint - look at the 'Total Size' in the output table"
                    },
                ],
            },
        ])

    def Q4(self):
        display_quiz([{
            "question":
                "What is the site-based Fst between population 0 and population 3 "
                "(to 4 decimal places)?",
            "type": "numeric",
            "precision": 4,
            "answers": [
                {
                    "type": "value",
                    "value": round(float(
                        self.mts.Fst([
                            self.mts.samples(population=0),
                            self.mts.samples(population=3)
                        ])
                    ), 4),
                    "correct": True,
                    "feedback":
                        "Correct"
                },
                {
                    "type": "range",
                    "range": [ -100000000, 1000000], 
                    "correct": False,
                    "feedback":
                        "Try again"
                },
            ]
        }])

    def Q5(self):
        display_quiz([{
            "question":
                "What is the mean branch-length Fst between samples from pop_0 and pop_1"
                " (to 4 decimal places)?",
            "type": "numeric",
            "precision": 4,
            "answers": [
                {
                    "type": "value",
                    "value": round(self.Fst_0_1_mean_100_reps, 4),
                    "correct": True,
                    "feedback":
                        "Correct"
                },
                {
                    "type": "range",
                    "range": [ -100000000, 1000000], 
                    "correct": False,
                    "feedback":
                        "Try again"
                },
            ]
        }])

    def Q6(self):
        display_quiz([
            {
                "question":
                    "The first site was used for inference; what is its inference_type"
                    " (as described in its metadata)?",
                "type": "multiple_choice",
                "answers": [
                    {"answer": "parsimony", "correct": False, "feedback": "Try again"},
                    {"answer": "full", "correct": True, "feedback": "Correct"},
                    {"answer": "fuzzy", "correct": False, "feedback": "Try again"},
                    {"answer": "potato", "correct": False, "feedback": "Really? You're joking, right?"},
                ]
            },
        ])

    def Q6bonus(self):
        display_quiz([
            {
                "question":
                    "In plot c, inferred using the default mismatch ratio of 1, how many"
                    " sites between 78 and 110 kb have been wrongly inferred to have"
                    " multiple mutations?",
                "type": "numeric",
                "precision": 0,
                "answers": [
                    {
                        "type": "value",
                        "value": 2,
                        "correct": True,
                        "feedback":
                            "Correct"
                    },
                    {
                        "type": "range",
                        "range": [ -100000000, 1000000], 
                        "correct": False,
                        "feedback":
                            "Try again"
                    },
                ]
            },
        ])


    def Q7(self):
        display_quiz([{
            "question":
                "What is the site-based genetic diversity in both the original and "
                "inferred tree sequences (to 5 decimal places)?",
            "type": "numeric",
            "precision": 5,
            "answers": [
                {
                    "type": "value",
                    "value": round(float(self.mts.diversity()), 5),
                    "correct": True,
                    "feedback":
                        "Correct: since this is a site-base measure, the same value will"
                        " be obtained regardless of how accurately the genealogy has"
                        " been inferred."
                },
                {
                    "type": "range",
                    "range": [ -100000000, 1000000], 
                    "correct": False,
                    "feedback":
                        "Try again"
                },
            ]
        }])

def setup_workbook1():
    return Workbook1()

def setup_workbook2():
    return Workbook2()
        

