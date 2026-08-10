import React from 'react';

import {selectModelNames, selectTemplate} from 'etna-js/selectors/magma';
import {useReduxState} from 'etna-js/hooks/useReduxState';
import ModelLink, {Arrowhead} from './model_link';
import ModelNode from './model_node';
import Layout from './tree_layout';
import Grid from '@material-ui/core/Grid';

const ModelMapGraphic = ({
  selected_models,
  handler,
  showLinks = false,
  width = 600,
  height = 600,
  disabled_models
}) => {
  let {model_names, templates} = useReduxState((state) => {
    let mod_names = selectModelNames(state);
    return {
      model_names: mod_names,
      templates: mod_names.map((model_name) =>
        selectTemplate(state, model_name)
      )
    };
  });

  let layout = new Layout(templates, width, height);

  return (
    <Grid style={{ width, height, position: 'relative' }}>
      <svg id='map' width={width} height={height}>
        <defs>
          <Arrowhead id='line_arrow' fill='#aca'/>
          <Arrowhead id='path_arrow' fill='goldenrod' fillOpacity='0.1'/>
        </defs>
        {showLinks && model_names.map((model_name) => {
          // for each link attribute, draw a link line.
          let node = layout.nodes[model_name];
          if (!node) return null;
          const { attributes } = node.model;

        return Object.entries(attributes).map( ([att_name, att]) => {
          if (att.attribute_type != 'link') return null;
        console.log({node, p: layout.nodes[att.link_model_name]});

          return <ModelLink
            key={`${model_name}.${att_name}`}
            center={ node.center}
            size={node.size}
            type='link'
            parent={ layout.nodes[att.link_model_name].center }
        />
        }).filter(_=>_);
        }).filter(_=>_).flat()}
        {model_names.map((model_name) => {
          let node = layout.nodes[model_name];
          return (
            <ModelLink
              key={model_name}
              center={node.center}
              parent={
                node.model.parent && layout.nodes[node.model.parent]
                  ? layout.nodes[node.model.parent].center
                  : null
              }
              size={node.size}
            />
          );
        })}

      </svg>
      {model_names.map((model_name) => {
        let node = layout.nodes[model_name];

        return (
          <ModelNode
            key={model_name}
            center={node.center}
            size={node.size}
            numSiblings={ layout.grid[node.depth] ? layout.grid[node.depth].length : 0 }
            selected={
              selected_models ? selected_models.includes(model_name) : false
            }
            handler={handler}
            model_name={model_name}
            disabled={
              disabled_models ? disabled_models.includes(model_name) : false
            }
          />
        );
      })}
    </Grid>
  );
};

export default ModelMapGraphic;
