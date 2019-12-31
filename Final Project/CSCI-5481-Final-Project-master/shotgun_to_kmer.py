from threading import Thread, Lock
from collections import defaultdict

THREADS = 8

### Take in shotgun data and output a list of k-mers
class ShotgunToKmer:

    def __init__(self, shotguns, kmer_length):
        self.kmer_length = kmer_length
        self.shotguns = shotguns
        self.shotgun_length = len(shotguns[0])  # assumes all shotguns are same length

        self.kmers = defaultdict(int)
        self.result_lock = Lock()

        threads = []
        num_shotguns_per_thread = len(shotguns) / THREADS
        # splits up shotgun kmer finding using threads
        for i in xrange(THREADS):
            start_shotgun_index = i*num_shotguns_per_thread
            if i == (THREADS - 1):  # add the remaining shotguns to last thread
                end_shotgun_index = len(shotguns)
            else:
                end_shotgun_index = (i+1)*num_shotguns_per_thread

            thread = Thread(target=self.find_kmers, args=(start_shotgun_index, end_shotgun_index))
            threads.append(thread)
            thread.start()

        for thread in threads:
            thread.join()

    def find_kmers(self, start_shotgun_index, end_shotgun_index):
        local_kmers = defaultdict(int)
        for shotgun_index in xrange(start_shotgun_index, end_shotgun_index):
            for start_index in xrange((self.shotgun_length-self.kmer_length)+1):
                end_index = start_index+self.kmer_length
                kmer = self.shotguns[shotgun_index][start_index:end_index]
                local_kmers[kmer] += 1

            if (shotgun_index - start_shotgun_index) % 10000:
                self.merge_kmers(local_kmers)
                local_kmers = defaultdict(int)

        self.merge_kmers(local_kmers)

    def merge_kmers(self, local_kmers):
        self.result_lock.acquire()

        for kmer, count in local_kmers.iteritems():
            self.kmers[kmer] = 1

        self.result_lock.release()

if __name__ == "__main__":
    # test case
    shotgun_inst = ShotgunToKmer(["zxcvmnasdf", "asdfqwerty", "tyuiopoiuy", "quickbrown"], kmer_length=3)
    print shotgun_inst.kmers
