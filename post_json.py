# importing the requests library 
import requests 
import json
import pprint

# defining the api-endpoint  
server_url = "http://127.0.0.1:5000/json/"
  
# your source code here 
post_dict_test1 = {"func" : "algebra.add", "args" : [1,3]}
post_dict_tether1 = {"func" : "tether",
                     "requests" : [{"func" : "algebra.add",
                                    "args" : [1,6]},
                                   {"func" : "algebra.sub",
                                    "args" : [3]},
                                   {"func" : "algebra.log",
                                    "args" : [2]}]}

post_json = json.dumps(post_dict_tether1)

# data to be sent to api 
data = post_json
  
# sending post request and saving response as response object 
r = requests.post(url = server_url, data = data) 

# extracting response text  
pprint.pprint(json.loads(r.text))
print("return type: ", type(json.loads(r.text)))