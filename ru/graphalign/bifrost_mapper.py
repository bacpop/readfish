import query_cpp

class Mapper:
    def __init__(self, index):
        self.index = index
        if self.index:
            self.mapper = query_cpp.Graph()
            self.mapper.read(index)
            self.initialised = True
        else:
            self.mapper = None
            self.initialised = False

    def map_read(self, seq, gap):
        return self.mapper.query(seq, gap)

    def map_reads(self, calls, gap):
        for read_id, seq in calls:
            yield read_id, list(self.mapper.query(seq, gap))

    def map_reads_2(self, calls, gap):
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
            yield read_info, read_id, seq_len, self.mapper.query(seq, gap)
