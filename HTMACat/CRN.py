import requests

def runCRN_net(configfname='ReactInfo'):
    with open(configfname) as f:
        config = f.read()
    post_dict = {'submit':'api', 'ReactInfo':config}
    url = 'http://www.catalysthub.net/CRN-api.php'
    r = requests.post(url, data=post_dict, verify=False)
    if r.status_code != 200:
        print('ERROR: status not 200 but', r.status_code)
        return
    else:
        print('Success...')
        return r.text
