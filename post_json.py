# importing the requests library 
import requests 
import json
import pprint

# defining the api-endpoint  
server_url = "http://127.0.0.1:5000/json/"
  
# your source code here 
post_dict_str = {"func" : "add", "args" : ["__help__"]}

post_dict_mixed_tether = {"func" : "tether",
             "args" : [
                       {"func" : "pipe1",
                            "args" : [[1,1,1],[1,2,3]]},
                       {"func" : "pipe2"},
                       {"func" : "tether",
                        "requests" : [{"func" : "pipe3"}]}
                       ]
                        
                      
            }

post_json = json.dumps(post_dict_str)

# data to be sent to api 
data = post_json
  
# sending post request and saving response as response object 
r = requests.post(url = server_url, data = data) 

# extracting response text  
pprint.pprint(json.loads(r.text))
print("return type: ", type(json.loads(r.text)))