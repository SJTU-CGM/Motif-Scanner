import math
from collections import OrderedDict


class MotifScanner:
    def __init__(self, sequence):
        self.sequence = sequence

    def scan(self, p_value, background, psdm):
        def windows(sequence, window_size):
            for i in range(0, len(sequence) - window_size + 1):
                yield i, sequence[i:i + window_size]

        def piece_is_valid(piece):
            for x in piece:
                if x not in alphabet:
                    return False
            return True

        alphabet = background.keys()
        pssm = psdm_to_pssm(psdm, background)
        score_distribution = get_score_distribution(background, pssm)
        score_threshold = get_score_threshold(score_distribution, p_value)
        for position, piece in windows(self.sequence, len(psdm)):
            if piece_is_valid(piece):
                piece_score = get_score(pssm, piece)
                if piece_score > score_threshold:
                    yield get_p_value(score_distribution, piece_score), position, piece


def psdm_to_pssm(psdm, background):
    # Position-Specific Probability Matrix --> Position-Specific Score Matrix
    pssm = []
    for distribution in psdm:
        code_to_score = {}
        for code, probability in distribution.items():
            if probability == 0:
                code_to_score[code] = float("-inf")
            else:
                code_to_score[code] = math.log2(probability / background[code])
        pssm.append(code_to_score)
    return pssm


def get_score_distribution(background, pssm):
    def update(d, k: float, p):
        d.setdefault(k, 0)
        d[k] += p
        return

    score_distribution = {0: 1}
    for score_table in pssm:
        dist = {}
        for base, p_base in background.items():
            for score, p_score in score_distribution.items():
                update(dist,
                       score + score_table[base],
                       p_score * p_base)
        score_distribution = dist
    return OrderedDict(sorted(score_distribution.items(), key=lambda x: x[0], reverse=True))


def get_score_threshold(distribution, p_value):
    last_score = None
    cumulative_probability = 0
    for score, probability in distribution.items():
        if cumulative_probability + probability > p_value:
            break
        else:
            last_score = score
            cumulative_probability += probability
    if last_score is None:
        raise Exception("You cannot search the instance of this pattern, as it is too similar to the background.")
    else:
        return last_score


def get_score(pssm, piece):
    score = 0
    for base_to_score, base in zip(pssm, piece):
        score += base_to_score[base]
    return score


def get_p_value(score_distribution: OrderedDict, piece_score):
    p_value = 0
    for score, probability in score_distribution.items():
        p_value += probability
        if piece_score == score:
            break
    return p_value
