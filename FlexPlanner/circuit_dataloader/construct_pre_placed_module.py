from collections import defaultdict

def construct_preplaced_modules(num_preplaced_module:int, blk_wh_dict:dict, outline_height:int) -> dict:
    """
    Add 'preplaced', 'x', 'y' to each block in blk_wh_dict.
    """
    for block_name in blk_wh_dict:
        blk_wh_dict[block_name]['preplaced'] = False

    blk_each_layer = defaultdict(list)
    for block_name in blk_wh_dict.keys():
        block_wh = blk_wh_dict[block_name]
        blk_each_layer[block_wh['z']].append({
            "name": block_name,
            "area": block_wh['w'] * block_wh['h'],
        })
    for i in range(len(blk_each_layer)):
        blk_each_layer[i] = sorted(blk_each_layer[i], key=lambda x: x["area"])
        x = 0
        for j,block in enumerate(blk_each_layer[i][0:num_preplaced_module]):
            block_name = block["name"]
            x += blk_wh_dict[block_name]['w']
            y = outline_height/2 if j >= num_preplaced_module//2 else 0

            blk_wh_dict[block_name]['preplaced'] = True
            blk_wh_dict[block_name]['x'] = x
            blk_wh_dict[block_name]['y'] = y

    return blk_wh_dict

