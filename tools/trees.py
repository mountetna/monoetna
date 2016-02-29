def dict_node(self, node, labels, default_name=None):
    if default_name:
        name = default_name
    elif node.id < len(labels):
        name = labels[node.id]
    else:
        name = 'internal_node'

    children = []

    if node.left:
      children.append( dict_node(node.left, labels) )
    if node.right:
      self.children.append( dict_node(node.right, labels) )

    return {
        'children': children,
        'name': name,
        'branch_length': node.dist
    }
