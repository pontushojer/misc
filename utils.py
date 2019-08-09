
def hamming_distance(s1, s2):
    # From https://pythonadventures.wordpress.com/2010/10/19/hamming-distance/
    assert len(s1) == len(s2)
    return sum(ch1 != ch2 for ch1, ch2 in zip(s1, s2))
