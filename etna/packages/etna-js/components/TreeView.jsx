import React, {useState, useEffect} from 'react';
import {CheckBox} from "./CheckBox";
import CollapsibleArrow from "./CollapsibleArrow";
import {trace} from "../utils/asserts";

require('./TreeView.css')

export default function TreeView({
  options = [],
  selected: selectedProp,
  disabled = {},
  collapsible = false,
  onSelectionsChange = (newState) => null,
  ItemLabelView = TreeViewDefaultItemLabelView,
  flowHorizontally = false,
} = {}) {
  const [selectedState, setSelected] = useState(() => initializeState(options));

  useEffect(() => {
    if (selectedProp == null) setSelected(initializeState(options));
  }, [options, selectedProp])

  const selected = selectedProp == null ? selectedState : selectedProp;

  return <div className={ flowHorizontally ? 'etna-tree-view vert' : 'etna-tree-view'}>
    {Nodes(options, selected, disabled)}
  </div>;

  function Nodes(options, selected, disabled={}, path = []) {
    const inner = (node, options, selected = {}, disabled = {}, path =[]) => WrapNode(
      <ItemLabelView node={node} path={path}
                     onChange={() => select(node, path, options, selected, disabled)}
                     disabled={isDisabled(options, disabled)}
                     checked={isSelected(options, selected)}/>,
      Nodes(options, selected, disabled, path.concat([node])));

    return options.map(([node, childOptions = []]) => <div className='child' key={node}>
      {inner(node, childOptions, selected[node], disabled[node], path)}
    </div>)
  }

  function WrapNode(label, children) {
    return collapsible ? <CollapsibleArrow label={label}>{children}</CollapsibleArrow> : <React.Fragment>
      <div>{label}</div>
      {children}
    </React.Fragment>
  }

  function isSelected(options, selected) {
    if (options.length === 0) return selected === true;
    return !options.find(([optionNode, childOptions = []]) => !isSelected(childOptions, selected[optionNode]))
  }

  function isDisabled(options, disabled) {
    if (options.length === 0) return disabled === true;
    return !options.find(([optionNode, childOptions = []]) => !isDisabled(childOptions, disabled[optionNode]))
  }


  function select(node, path, nodeOptions, nodeSelected, nodeDisabled) {
    const newState = { ...selected };
    const newSelected = path.reduce((s, p) => (s[p] = { ...s[p] }), newState);
    const originalSelected = path.reduce((s, p) => s[p], selected);

    if (nodeDisabled !== true) {
      newSelected[node] = initializeState(nodeOptions, node, !isSelected(nodeOptions, nodeSelected), nodeDisabled);

      function restoreDisabled(disabled = {}, original = {}, selected = {}) {
        Object.keys(disabled).forEach(k => {
          if (disabled[k] === true) selected[k] = original[k]
          if (disabled[k] != null) restoreDisabled(disabled[k], original[k], selected[k]);
        });
      }

      console.log(nodeDisabled, node, originalSelected, newSelected);
      restoreDisabled(nodeDisabled, originalSelected[node], newSelected[node]);
    }

    setSelected(newState);
    onSelectionsChange(newState);
  }
}

export function TreeViewDefaultItemLabelView({ node, onChange, checked, path, disabled }) {
  return <label>
    <CheckBox data-node={node} data-path={path} onChange={onChange} checked={checked} disabled={disabled} />
    {node}
  </label>
}

function initializeState(options, parent = undefined, leafValue = true) {
  if (options.length === 0 && parent !== undefined) return leafValue;
  return options.reduce((o, [node, options = []]) => (o[node] = initializeState(options, o, leafValue), o), {});
}

export function getSelectedLeaves(selectedState) {
  return Object.keys(selectedState).reduce((r, k) => r.concat(
    selectedState[k] === true ? [k] : getSelectedLeaves(selectedState[k])), []);
}