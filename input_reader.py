import configparser

def read_input_file(input_file):
    config = configparser.ConfigParser()
    config.read(input_file)

    # Step 1: Read fixed fields
    ion_type_num = int(config["Ion Type Number"]["ion_type_num"])
    solute_type = config["Solute Information"]["solute_type"]
    charge_rescale = config["Charge Rescale Factor"]["charge_rescale"]

    # Step 2: Read ion types dynamically
    ion_types = []
    for i in range(1, ion_type_num + 1):
        key = f"ion_type{i}"
        ion_types.append(config["Ion Type Information"][key])
        
    solute_type_raw = config["Solute Information"].get("solute_type", "").strip()
    solute_type = solute_type_raw if solute_type_raw else None


    parameters = {
        "ion_type_num": ion_type_num,
        "ion_types": ion_types,
        "solute_type": solute_type,
        "charge_rescale": charge_rescale
    }

    return parameters