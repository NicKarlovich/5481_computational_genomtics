import string, random, copy
from threading import Thread, Lock

THREADS = 8

### take in a sequence and output simulated shotgun sequences
class SimulatedShotgun:

    def __init__(self, sequence, shotgun_length, coverage, error_percent=0):
        self.sequence = sequence
        self.shotgun_length = shotgun_length

        self.shotguns = []
        self.result_lock = Lock()

        self.seq_length = len(sequence)

        # determine the total number of shotgun strings to return
        self.num_shotguns = (coverage * self.seq_length) / shotgun_length

        threads = []
        for i in xrange(THREADS):
            shotguns_to_make = self.num_shotguns / THREADS
            if i == (THREADS - 1):  # add the remaining shotguns to last thread
                shotguns_to_make += self.num_shotguns % THREADS

            thread = Thread(target=self.get_shotguns, args=[shotguns_to_make])
            threads.append(thread)
            thread.start()

        for thread in threads:
            thread.join()

        self.simulated_errors(error_percent)

    def simulated_errors(self, error_percent):
        if error_percent:
            num_errors = int(self.num_shotguns*self.shotgun_length * (error_percent / 100.0))
            dna_bases = ["a", "g", "t", "c"]
            for _ in xrange(num_errors):
                rand_shotgun = random.randint(0, self.num_shotguns-1)
                rand_index = random.randint(0, self.shotgun_length-1)
                rand_selection = self.shotguns[rand_shotgun][rand_index]

                adj_bases = copy.deepcopy(dna_bases)
                adj_bases.remove(rand_selection.lower())

                rand_base_index = random.randint(0,len(adj_bases)-1)
                self.shotguns[rand_shotgun] = self.shotguns[rand_shotgun][:rand_index] + adj_bases[rand_base_index] + self.shotguns[rand_shotgun][rand_index+1:]

    def get_shotguns(self, shotguns_to_make):
        local_shotguns = []
        max_index = self.seq_length - self.shotgun_length

        for _ in xrange(shotguns_to_make):
            # pick a random start position
            start_index = random.randint(0, max_index)
            shotgun = self.sequence[start_index: start_index+self.shotgun_length]
            local_shotguns.append(shotgun)

        self.result_lock.acquire()
        self.shotguns += local_shotguns
        self.result_lock.release()


if __name__ == "__main__":
    # just a test with the alphabet
    test_str = ''.join(["AGTCGCATAGGCCATT"]*10)
    shotgun_inst = SimulatedShotgun(test_str, shotgun_length=10, coverage=10, error_percent=3)
    print shotgun_inst.shotguns
