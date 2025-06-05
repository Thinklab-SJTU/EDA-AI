import torch
# from evaluate.preprocess.utils import output


class SteinerExtractor():
    def __init__(self):
        pass

    def GetAllPoints(self,points_maps:torch.Tensor)->torch.Tensor:

        output_tensor = torch.full(points_maps.size(),-1.0)
        output_tensor[:,1,:,:] = points_maps[:,1,:,:]
        output_tensor[:,2,:,:] = points_maps[:,2,:,:]
        mask = points_maps[:, 0, :, :] >= 0.8
        output_tensor[:, 0, :, :][mask] = 1.0
        # output_tensor[:, 0, :, :] = torch.where(mask, torch.tensor(1.0, dtype=points_maps.dtype), torch.tensor(-1.0, dtype=points_maps.dtype))
        return output_tensor

    def GetPredictedPoints(self,points_maps,conditions,meta_value = torch.tensor(1.0))->torch.Tensor:
        mask = conditions[:,0,:,:] == meta_value
        output_tensor = torch.full(points_maps.size(),-1.0)
        output_tensor[:, 1, :, :] = points_maps[:, 1, :, :]
        output_tensor[:, 2, :, :] = points_maps[:, 2, :, :]
        output_tensor[:, 0, :, :][mask] = -1.0
        return output_tensor

    def GetStenierPoints(self,conditions,predictedpoints)->torch.Tensor:
        pass


if __name__ == '__main__':
    extrator = SteinerExtractor()
    input_tensor = torch.randn((3,3,4,4))
    output_tensor = extrator.GetAllPoints(input_tensor)
    print(f'input_tensor:{input_tensor}')
    print(f'output_tensor:{output_tensor}')
