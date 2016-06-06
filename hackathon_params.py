import json


def hackParams(N, stp):
    
    params = {'N' : N,  'steps' : stp}

    with open('params.json', 'w') as fp:
        json.dump(params, fp)
