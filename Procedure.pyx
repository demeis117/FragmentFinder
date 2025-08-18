
# Procedure_sorted_uid_match_fixed_returning_total.pyx
# cython: language_level=3

import numpy as np
cimport numpy as cnp
from ahocorasick import Automaton

cpdef dict original_strings(str db_file):
    cdef dict sequences = {}
    cdef str header = ""
    cdef list seq_parts = []
    cdef str line
    with open(db_file, 'r') as fh:
        for line in fh:
            if line.startswith('>'):
                if header:
                    sequences[header] = "".join(seq_parts)
                header = line[1:].strip().split()[0]
                seq_parts = []
            else:
                seq_parts.append(line.strip())
        if header:
            sequences[header] = "".join(seq_parts)
    return sequences


cpdef separate_seqs(str db_file, str fasta_file, str out_file, int k):
    cdef dict ref_seqs = original_strings(db_file)
    cdef object A = Automaton()
    cdef str seq, line, hdr, read_seq
    cdef int pos, total_in = 0, total_out = 0

    for seq in ref_seqs.values():
        if len(seq) < k:
            continue
        for pos in range(len(seq) - k + 1):
            A.add_word(seq[pos:pos+k], True)
    A.make_automaton()

    hdr = ""
    read_seq = ""

    with open(fasta_file, 'r') as inp, open(out_file, 'w') as outp:
        for line in inp:
            line = line.strip()
            if line.startswith('>'):
                if hdr and read_seq:
                    total_in += 1
                    try:
                        if next(A.iter(read_seq), None) is not None:
                            outp.write(hdr + "\n" + read_seq + "\n")
                            total_out += 1
                    except:
                        pass
                hdr = line
                read_seq = ""
            else:
                read_seq += line.strip()
        if hdr and read_seq:
            total_in += 1
            try:
                if next(A.iter(read_seq), None) is not None:
                    outp.write(hdr + "\n" + read_seq + "\n")
                    total_out += 1
            except:
                pass

    print(f"[separate_seqs] Total reads processed: {total_in}")
    print(f"[separate_seqs] Reads matched to k-mers: {total_out}")
    return total_in


cpdef tuple seq_aligner(str db_file, str reads_file, int k):
    cdef dict ref_seqs = original_strings(db_file)
    cdef list headers = list(ref_seqs.keys())
    cdef dict raw_results_dict = {}
    cdef dict metadata_dict = {}
    cdef object A = Automaton()
    cdef int idx, pos, total_reads = 0, end_index, peak_val
    cdef str uid, ref_seq, line, read_seq
    cdef tuple payload
    cdef double norm_val
    cdef cnp.ndarray[cnp.int64_t, ndim=1] vec

    for idx, uid in enumerate(headers):
        ref_seq = ref_seqs[uid]
        vec = np.zeros(len(ref_seq), dtype=np.int64)
        raw_results_dict[uid] = vec
        if len(ref_seq) >= k:
            for pos in range(len(ref_seq) - k + 1):
                A.add_word(ref_seq[pos:pos+k], (idx, pos))
    A.make_automaton()

    read_seq = ""

    with open(reads_file, 'r') as rf:
        for line in rf:
            if line.startswith('>'):
                if read_seq:
                    total_reads += 1
                    for end_index, payload in A.iter(read_seq):
                        idx, pos = payload
                        uid = headers[idx]
                        raw_results_dict[uid][pos] += 1
                read_seq = ""
            else:
                read_seq += line.strip()

        if read_seq:
            total_reads += 1
            for end_index, payload in A.iter(read_seq):
                idx, pos = payload
                uid = headers[idx]
                raw_results_dict[uid][pos] += 1

    for uid in headers:
        vec = raw_results_dict[uid]
        peak_val = int(vec.max()) if vec.size else 0
        norm_val = peak_val * 1e6 / total_reads if total_reads > 0 else 0.0
        metadata_dict[uid] = norm_val

    return raw_results_dict, metadata_dict
