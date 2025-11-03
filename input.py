import json

def read_json(json_name):
    try:
        with open(f"Examples/{json_name}.json","r") as f:
            data = json.load(f)
            return data
    except Exception as e :
        print(f"Couldn't read the json file because : {e}")
