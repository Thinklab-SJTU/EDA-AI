from collections import OrderedDict

def assign_layer(blk_wh_dict:dict, num_layer:int) -> OrderedDict:
    """
    Assign layer for each block. 
    Aiming to make the total area of each layer as close as possible.
    Add 'z' key to each block in blk_wh_dict, which represents the layer index.
    """
    # after sorting, blk_items is a list of items, the first is the block name, the second is the block info
    blk_items = sorted(blk_wh_dict.items(), key=lambda x: x[1]['w'] * x[1]['h'], reverse=True)

    item_each_layer = [ [] for _ in range(num_layer) ]
    area_sum_each_layer = [ 0 for _ in range(num_layer) ]
    
    for i in range(len(blk_items)):
        # add current item to the layer with the smallest area sum
        min_idx = area_sum_each_layer.index(min(area_sum_each_layer))
        blk_items[i][1]['z'] = min_idx
        item_each_layer[min_idx].append(blk_items[i])
        area_sum_each_layer[min_idx] += blk_items[i][1]['w'] * blk_items[i][1]['h']

    return OrderedDict(sum(item_each_layer, []))
