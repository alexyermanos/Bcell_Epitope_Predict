from torch_geometric.nn import GATConv
import torch.nn.functional as F
import torch

class GAT(torch.nn.Module):
    def __init__(self, in_dim, hid_dim, out_dim, num_head, out_head):
        super(GAT, self).__init__()
        
        self.in_dim = in_dim
        self.hid_dim = hid_dim
        self.out_dim = out_dim
        
        self.in_head = num_head
        self.out_head = out_head

        self.conv_in = GATConv(self.in_dim, self.hid_dim, heads=self.in_head) 
        self.conv_mid = GATConv(self.hid_dim*self.in_head, self.hid_dim, heads=self.in_head)
        self.conv_out = GATConv(self.hid_dim*self.in_head, self.out_dim, concat=False, heads=self.out_head)
        self.sigmoid = torch.nn.Sigmoid()
        
        
    def forward(self, data):
        
        x, edge_index = data.x, data.edge_index
        
        x = self.conv_in(x, edge_index) 
        initial_x = x
        x = F.elu(x)
        
        x = self.conv_mid(x, edge_index) + initial_x
        x = F.elu(x)
        
        x = self.conv_mid(x, edge_index) + initial_x
        x = F.elu(x)
        
        x = self.conv_mid(x, edge_index) + initial_x
        x = F.elu(x)
        
        x = self.conv_mid(x, edge_index) + initial_x
        x = F.elu(x)
        
        x = self.conv_mid(x, edge_index) + initial_x
        x = F.elu(x)
        
        x = self.conv_mid(x, edge_index) + initial_x
        x = F.elu(x)
        
        x = self.conv_out(x, edge_index) 
        x = self.sigmoid(x)
        
        return x