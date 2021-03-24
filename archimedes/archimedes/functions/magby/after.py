import re

# We should probably improve this function for cases when the return is more complicated
def query_extract(data, target):
    out=data['answer']
    if 'format' in data.keys():
        ind=1
        out=[]
        for i in data['answer']:
            out = out + [i[ind]]
    return out
