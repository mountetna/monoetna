def add_node(node,tree_dict,labels):
    '''
    Create a nested dictionary from the ClusterNode's returned by SciPy
    '''
    # create the new node and append it to its parent's children
    if node.id < len(labels):
        name = labels[node.id]
    else:
        name= 'internal_node'
    new_node = {'name': name , 'children':[], 'branch_length': node.dist}
    tree_dict['children'].append( new_node )
    # Recursively add the current node's children
    if node.left:
        add_node(node.left, new_node, labels)
    if node.right:
        add_node(node.right, new_node, labels)
    return tree_dict
    
def main():
    
    print 'This is tree_methods.py'
    

if __name__ == "__main__":
    main()
    