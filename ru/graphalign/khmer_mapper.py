import khmer
from khmer import khmer_args
from .long_align import align_long

class Mapper:
    def __init__(self, index):
        self.index = index
        if self.index:
            self.graph = khmer.Countgraph.load(index)
            self.mapper = khmer.ReadAligner(self.graph, 1)
            self.initialised = True
        else:
            self.mapper = None
            self.graph = None
            self.initialised = False

    def map_read(self, seq):
        return align_long(self.index, self.mapper, seq)

    def map_reads(self, calls):
        for read_id, seq in calls:
            yield read_id, list(align_long(self.index, self.mapper, seq))

    def map_reads_2(self, calls):
        """Align reads against a reference
        Parameters
        ----------
        calls : iterable [tuple,  str, str, int, str]
            An iterable of called reads from PerpetualCaller.basecall_minknow
        Yields
        ------
        read_info : tuple
            Tuple of read info (channel, read_number)
        read_id : str
        sequence : str
        sequence_length : int
        mapping_results : list
        """
        for read_info, read_id, seq, seq_len, quality in calls:
            yield read_info, read_id, seq_len, align_long(self.graph, self.mapper, seq)
