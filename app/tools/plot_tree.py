import os
from dotenv import load_dotenv
from Bio import Phylo
import pylab
import networkx

load_dotenv()
DATA_PATH = os.getenv('DATA_PATH')

tree = Phylo.read(os.path.join(DATA_PATH, 'sequence_tree.dnd'), 'newick')
# tree.ladderize()
net = Phylo.to_networkx(tree)
networkx.draw(net)
pylab.show()
