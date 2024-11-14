from matplotlib import pyplot as plt
from fp_env import FPInfo
from fp_env import Block
import torch
import math
import numpy as np
import pandas as pd
from collections import defaultdict

def save_final_floorplan(path:str, fp_info:FPInfo, enable_text:bool=1) -> None:
    """
    draw and save the final floorplanning result.
    ALso save a .csv file with the same name as the image file.
    """
    num_layer = fp_info.num_layer
    label_dict = {
        "Preplaced": [0] * num_layer,
        "Movable": [0] * num_layer,
    }

    df = pd.DataFrame()
    name2alignment_group = fp_info.name2alignment_group

    fig, axes = plt.subplots(1, num_layer, figsize=(10, 6))
    if isinstance(axes, plt.Axes):
        axes = [axes]
    for block in fp_info.block_info:
        if block.virtual:
            continue
        
        if block.placed:
            facecolor = "gray" if block.preplaced else fp_info.name2alignment_group_color[block.name]
            label = "Preplaced" if block.preplaced else "Movable"
            label_dict[label][block.grid_z] += 1
            label = label if label_dict[label][block.grid_z] == 1 else None

            axes[block.grid_z].add_patch(plt.Rectangle((block.grid_x, block.grid_y), block.grid_w, block.grid_h, fill=True, facecolor=facecolor, edgecolor="k", alpha=0.3, label=label))
            # text = str(block.idx)
            text = block.name
            if block.name in name2alignment_group:
                text += " (A {})".format(name2alignment_group[block.name])
            if enable_text:
                axes[block.grid_z].text(block.grid_x + block.grid_w/2, block.grid_y + block.grid_h/2, text, ha='center', va='center')
            df = pd.concat([df, pd.DataFrame([{
                "name": block.name,
                "x": block.grid_x,
                "y": block.grid_y,
                "z": block.grid_z,
                "w": block.grid_w,
                "h": block.grid_h,
                "type": "block",
                "preplaced": block.preplaced,
            }])], ignore_index=True)
    
    # draw terminal
    terminal_coordinates_x = [[] for _ in range(num_layer)]
    terminal_coordinates_y = [[] for _ in range(num_layer)]
    for terminal in fp_info.terminal_info:
        terminal_coordinates_x[terminal.grid_z].append(terminal.grid_x)
        terminal_coordinates_y[terminal.grid_z].append(terminal.grid_y)
        df = pd.concat([df, pd.DataFrame([{
            "name": terminal.name,
            "x": terminal.grid_x,
            "y": terminal.grid_y,
            "z": terminal.grid_z,
            "type": "terminal",
        }])], ignore_index=True)
    for i in range(num_layer):
        if len(terminal_coordinates_x[i]) > 0:
            axes[i].scatter(terminal_coordinates_x[i], terminal_coordinates_y[i], color="red", label="Terminal", s=5)


    # draw outline
    chip_w = fp_info.x_grid_num
    chip_h = fp_info.y_grid_num
    outline = max(chip_w, chip_h)

    for i in range(num_layer):
        axes[i].set_aspect('equal')
        axes[i].set_xlim([-0.1*outline, outline * 1.1])
        axes[i].set_ylim([-0.1*outline, outline * 1.1])

        axes[i].plot([chip_w, chip_w], [0, chip_h], color='k')
        axes[i].plot([0, chip_w], [chip_h, chip_h], color='k')
        axes[i].plot([0, 0], [0, chip_h], color='k')
        axes[i].plot([0, chip_w], [0, 0], color='k')

        # set title
        axes[i].set_title(f"Layer {i}")
        axes[i].legend(loc="upper right")


    # save
    hpwl, weighted_hpwl, original_hpwl = fp_info.calc_hpwl()
    title = "HPWL = {}, original HPWL = {}".format( hpwl, int(original_hpwl) )
    alignment = fp_info.calc_alignment_score()
    overlap = fp_info.get_overlap()
    title += "\nAlignment rate = {:.2f}, Overlap = {:.2f}".format(alignment, overlap)
    plt.suptitle(title)
    plt.savefig(path)
    plt.close()

    df = pd.concat([pd.DataFrame([{
        "hpwl": hpwl,
        "original_hpwl": original_hpwl,
        "overlap": overlap,
        "alignment": alignment,
    }]), df], ignore_index=True)
    df.to_csv(path.replace(".png", ".csv"), index=False)


def save_intermediate_floorplan(path:str, curr_block:Block, canvas:torch.Tensor,
              position_mask:torch.Tensor, wiremask:torch.Tensor, alignment_mask:torch.Tensor, binary_alignment_mask:torch.IntTensor, 
              fp_info:FPInfo, enable_text:bool=1,
              ) -> None:
    """draw and save all masks at current step."""
    num_layer = fp_info.num_layer
    label_dict = {
        "Preplaced": [0] * num_layer,
        "Movable": [0] * num_layer,
    }

    num_axes = num_layer # [num_layer] canvas with rectangle
    if canvas is not None:
        num_axes += num_layer
    if wiremask is not None:
        num_axes += 1
    if position_mask is not None:
        num_axes += 1
    if alignment_mask is not None:
        num_axes += 1
    if binary_alignment_mask is not None:
        num_axes += 1
    masks_for_available_mask = list(filter(lambda x: x is not None, [position_mask, binary_alignment_mask]))
    if len(masks_for_available_mask) > 0: # available mask
        num_axes += 1


    nrows = 4
    ncols = math.ceil(num_axes/nrows)
    size = 5
    fig, axes = plt.subplots(nrows, ncols, figsize=(size*ncols, size*nrows)) # W x H
    axes = axes.flatten()

    # draw block on canvas
    for block in fp_info.block_info:
        if block.virtual:
            continue
        if block.placed:
            facecolor = "gray" if block.preplaced else fp_info.name2alignment_group_color[block.name]
            label = "Preplaced" if block.preplaced else "Movable"
            label_dict[label][block.grid_z] += 1
            label = label if label_dict[label][block.grid_z] == 1 else None
            
            axes[block.grid_z].add_patch(plt.Rectangle((block.grid_x, block.grid_y), block.grid_w, block.grid_h, fill=True, facecolor=facecolor, edgecolor="k", alpha=0.3, label=label))
            if enable_text:
                axes[block.grid_z].text(block.grid_x + block.grid_w/2, block.grid_y + block.grid_h/2, str(block.idx),ha='center', va='center')
    
    # draw terminal
    terminal_coordinates_x = []
    terminal_coordinates_y = []
    for terminal in fp_info.terminal_info:
        if terminal.grid_z == curr_block.grid_z:
            terminal_coordinates_x.append(terminal.grid_x)
            terminal_coordinates_y.append(terminal.grid_y)
    if len(terminal_coordinates_x) > 0:
        axes[0].scatter(terminal_coordinates_x, terminal_coordinates_y, color="red", label="Terminal", s=5)
        

    # draw outline for each layer
    chip_w = fp_info.x_grid_num
    chip_h = fp_info.y_grid_num
    outline = max(chip_w, chip_h)

    for z in range(num_layer):
        axes[z].set_aspect('equal')
        axes[z].set_xlim([-0.1*outline, outline * 1.1])
        axes[z].set_ylim([-0.1*outline, outline * 1.1])

        axes[z].plot([chip_w, chip_w], [0, chip_h], color='k')
        axes[z].plot([0, chip_w], [chip_h, chip_h], color='k')
        axes[z].plot([0, 0], [0, chip_h], color='k')
        axes[z].plot([0, chip_w], [0, 0], color='k')

        # set title
        axes[z].set_title(f"Layer {z}")
        axes[z].legend(loc="upper right")


    # rotate position mask and draw with imshow
    axes_idx = num_layer
    
    # draw canvas of current layer
    if canvas is not None:
        for z in range(num_layer):
            canvas_z = canvas[z]
            canvas_z = torch.rot90(canvas_z, 1)
            x_max, y_max = canvas_z.shape
            plt.colorbar(axes[axes_idx].imshow(canvas_z, extent=[0, x_max, 0, y_max], interpolation='none'), ax=axes[axes_idx])
            axes[axes_idx].set_title(f"canvas, layer={z}")
            axes_idx += 1
            canvas_z = torch.rot90(canvas_z, -1)
    
    # draw wiremask
    if wiremask is not None:
        wiremask = torch.rot90(wiremask, 1)
        x_max, y_max = wiremask.shape
        plt.colorbar(axes[axes_idx].imshow(wiremask, extent=[0, x_max, 0, y_max], interpolation='none'), ax=axes[axes_idx])
        axes[axes_idx].set_title(f"wiremask, layer={curr_block.grid_z}\nfull_idx={curr_block.idx}, movable_idx={curr_block.movable_idx}")
        axes_idx += 1
        wiremask = torch.rot90(wiremask, -1)
    
    # draw position mask
    if position_mask is not None:
        position_mask = torch.rot90(position_mask, 1)
        x_max, y_max = position_mask.shape
        plt.colorbar(axes[axes_idx].imshow(position_mask, extent=[0, x_max, 0, y_max], interpolation='none'), ax=axes[axes_idx])
        axes[axes_idx].set_title(f"position_mask, layer={curr_block.grid_z}\nfull_idx={curr_block.idx}, movable_idx={curr_block.movable_idx}")
        axes_idx += 1
        position_mask = torch.rot90(position_mask, -1)
    
    # draw alignment mask
    if alignment_mask is not None:
        alignment_mask = torch.rot90(alignment_mask, 1)
        x_max, y_max = alignment_mask.shape
        plt.colorbar(axes[axes_idx].imshow(alignment_mask, extent=[0, x_max, 0, y_max], interpolation='none'), ax=axes[axes_idx])
        axes[axes_idx].set_title(f"alignment_mask, layer={curr_block.grid_z}\nfull_idx={curr_block.idx}, movable_idx={curr_block.movable_idx}\npartner_indices={curr_block.partner_indices}")
        axes_idx += 1
        alignment_mask = torch.rot90(alignment_mask, -1)
    
    # draw binary alignment mask
    if binary_alignment_mask is not None:
        binary_alignment_mask = torch.rot90(binary_alignment_mask, 1)
        x_max, y_max = binary_alignment_mask.shape
        plt.colorbar(axes[axes_idx].imshow(binary_alignment_mask, extent=[0, x_max, 0, y_max], interpolation='none'), ax=axes[axes_idx])
        axes[axes_idx].set_title(f"binary_alignment_mask, layer={curr_block.grid_z}\nfull_idx={curr_block.idx}, movable_idx={curr_block.movable_idx}\npartner_indices={curr_block.partner_indices}")
        axes_idx += 1
        binary_alignment_mask = torch.rot90(binary_alignment_mask, -1)

    
    # draw available mask
    if len(masks_for_available_mask) > 0:
        available_mask = torch.stack(masks_for_available_mask, dim=0).sum(dim=0)
        available_mask = torch.rot90(available_mask, 1)
        x_max, y_max = available_mask.shape
        plt.colorbar(axes[axes_idx].imshow(available_mask, extent=[0, x_max, 0, y_max], interpolation='none'), ax=axes[axes_idx])
        axes[axes_idx].set_title(f"available_mask, layer={curr_block.grid_z}\nfull_idx={curr_block.idx}, movable_idx={curr_block.movable_idx}")
        axes_idx += 1
        available_mask = torch.rot90(available_mask, -1)

    

    super_title = f"name = {curr_block.name}, layer = {curr_block.grid_z} \n full_idx = {curr_block.idx}, movable_idx = {curr_block.movable_idx} \n partner_indices = {curr_block.partner_indices}"
    plt.suptitle(super_title)

    # save
    plt.savefig(path)
    plt.close()



def draw_action_record(action_record:dict[str, dict[str, list]], save_path:str):
    """
    {
        "env_id": {
            "action": [action1, action2, ...],
        }
    }
    """
    num_env = len(action_record)
    assert num_env > 0
    num_kind_act = len(list(action_record.values())[0])
    assert num_kind_act > 0

    fig, axes = plt.subplots(num_env, num_kind_act, figsize=(5*num_kind_act, 5*num_env))
    if not isinstance(axes, np.ndarray):
        axes = np.array([axes])
    axes = axes.reshape(num_env, num_kind_act)
    
    for env_id, act_dict in action_record.items():
        for act_idx, (act_name, act) in enumerate(act_dict.items()):
            env_id = int(env_id)
            if act_name in ["layer"]:
                act = np.array(act).astype(np.int32)
                act = np.cumsum(act)
                # axes[env_id, act_idx].plot(act, 'o', label=act_name, markersize=0.5)
                axes[env_id, act_idx].plot(act, label=act_name)
                axes[env_id, act_idx].plot([0, len(act)], [0, len(act) // 2], label="cumsum", linestyle="--", color="red")
            else:
                axes[env_id, act_idx].plot(act, label=act_name)

            axes[env_id, act_idx].set_title(f"env_id={env_id}, act_name={act_name}")
            axes[env_id, act_idx].legend()
            # grid
            axes[env_id, act_idx].grid()
            # xlabel
            axes[env_id, act_idx].set_xlabel("step")
            # ylabel
            axes[env_id, act_idx].set_ylabel(act_name)

    plt.savefig(save_path)
    plt.close()