from constants import STRING_EXTRACTED_URL_PATH, HI_UNION_PATH, LIT_BM_PATH


def get_average_huri_score():
    '''
    Get the average score for HuRI interactions in the STRING dataset
    '''

    with open(STRING_EXTRACTED_URL_PATH, "r") as string_links:

        scoredict = {}
        for link in string_links.readlines()[1:]:
            p1, p2, score = link.strip().split(' ')
            p1 = p1.split('.')[1]
            p2 = p2.split('.')[1]
            score = int(score)
            scoredict['_'.join(sorted([p1, p2]))] = score

    huri_score_list = []
    with open(HI_UNION_PATH, "r") as huri_links, open(LIT_BM_PATH, "r") as lit_links:
        interactions = huri_links.readlines() + lit_links.readlines()
        for interaction in interactions:
            p1, p2 = interaction.strip().split()
            lookup_key = '_'.join(sorted([p1, p2]))
            lookup_score = scoredict.get(lookup_key, None)
            if lookup_score:
                huri_score_list.append(lookup_score)

    return huri_score_list

huri_score_list = get_average_huri_score()

