#!/usr/bin/env python3

import requests
import json
import argparse, textwrap


def collapse(jsondata):

    def goThrough(x):
        if isinstance(x, list):
            return goThroughList(x)
        elif isinstance(x, dict):
            return goThroughDict(x)
        else:
            return x

    def goThroughDict(dic):
        for key, value in dic.items():
            if isinstance(value, dict):
                dic[key] = goThroughDict(value)
            elif isinstance(value, list):
                dic[key] = goThroughList(value)
        return dic

    def goThroughList(lst):
        if(not any(not isinstance(y, (int, float)) for y in lst)):  # pure numeric list
            if len(lst) <= 2:
                return lst
            else:
                return '['+str(lst[0]) + ' ... '+str(lst[-1])+'] ('+str(len(lst))+')'
        else:
            return [goThrough(y) for y in lst]

    return goThroughDict(jsondata)


def callEtc(postdatafile, url, uploadfile=None):

    with open(postdatafile) as f:
        postdata = json.loads(f.read())

    # TODO! workaround until we figure put how to handle ssl certificate correctly
    import warnings
    warnings.filterwarnings('ignore', message='Unverified HTTPS request')

    if uploadfile is None:
        return requests.post(url,
                          data=json.dumps(postdata),
                          headers={'Content-Type': 'application/json'},
                          verify=False,
                          )
        return r
    else:
        return requests.post(url,
                          data={"data": json.dumps(postdata)},
                          files={"target": open(uploadfile, 'rb')},
                          verify=False)

def output(jsondata, args):

    if(args.collapse):
        jsondata = collapse(jsondata)

    if(args.outputfile):
        with open(args.outputfile, "w") as of:
            of.write(json.dumps(jsondata, indent=args.indent))
    else:
        print(json.dumps(jsondata, indent=args.indent))


def getEtcUrl(etcname):
    if('4most' in etcname.lower()):
        return '4Most/'
    elif('crires' in etcname.lower()):
        return 'Crires2/'
    # elif('harps' in etcname.lower()):  # not yet working
    #     return 'HarpsNirps/'
    else:
        print("error: no match for etcname: " + etcname)


def main():
    parser = argparse.ArgumentParser(
        description='Call an ETC with input parameters and optionally an uploaded spectrum.\n'
        + 'Print the resulting JSON on stdout or optionally a file.\n'
        + 'Examples: \n'
        + './etc_cli.py crires etc-form.json -o output1.json\n'
        + './etc_cli.py crires etc-form-uploading.json -u upload.dat -o output2.json',
        formatter_class=argparse.RawTextHelpFormatter
    )

    parser.add_argument('etcname',
                        help='Name of instrument ETC to call, e.g. 4most')

    parser.add_argument('postdatafile',
                        help='Name of JSON file with ETC input parameters,\nlike the ETC input form')

    parser.add_argument('-u,', '--upload', dest="uploadfile",
                        help='Name of file with spectrum to upload.\nSee https://etc.eso.org/observing/etc/doc/upload.html')

    parser.add_argument('-c', '--collapse', action='store_true',
                        help='collapse output JSON data arrays to a short indicative strings')

    parser.add_argument('-i', '--indent', type=int, nargs='?', const=4,
                        help='Format the output JSON with indentation (default 4)')

    parser.add_argument('-o', '--outputfile', dest="outputfile",
                        help='Send the output to file')

    args = parser.parse_args()

    # baseurl = 'http://localhost:8000/observing/etc/etcapi/'
    baseurl = 'https://etctest.hq.eso.org/observing/etc/etcapi/'
    # baseurl = 'https://etc.eso.org/observing/etc/etcapi/'

    etcName = getEtcUrl(args.etcname)

    url = baseurl + getEtcUrl(etcName)

    jsondata = callEtc(args.postdatafile, url, args.uploadfile).json()

    output(jsondata, args)


if __name__ == "__main__":
    main()
