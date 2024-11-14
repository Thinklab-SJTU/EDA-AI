class Terminal:
    def __init__(self, x:float, y:float, z:int, name:str) -> None:
        self.x = x
        self.y = y
        self.z = z
        self.name = name

        self.grid_z = z

        self.placed = True
        self.connected_nets = []
    
    def __repr__(self) -> str:
        return "Terminal(name={}, id={}, terminal_idx={}, grid_x={}, grid_y={}, grid_z={}, x={}, y={})".format(
            self.name,
            self.idx, self.terminal_idx, self.grid_x, self.grid_y, self.grid_z,
            self.x, self.y
        )


    def set_grid_xy(self, grid_width:float, grid_height:float, x_grid_num:int, y_grid_num:int):
        """set grid_x and grid_y based on x, y, grid_width and grid_height."""
        self.grid_x = round(min(x_grid_num - 1, self.x / grid_width))
        self.grid_y = round(min(y_grid_num - 1, self.y / grid_height))


    def set_idx(self, full_idx:int, terminal_idx:int):
        self.idx = full_idx
        self.terminal_idx = terminal_idx
