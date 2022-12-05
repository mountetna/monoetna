// Framework libraries.
import * as React from 'react';

import {MagmaName} from './magma_item';
import {DefaultName} from '../item_name';

import MarkdownAttribute from 'etna-js/components/attributes/markdown_attribute';
import AttributeViewer from 'etna-js/components/attributes/attribute_viewer';
import MarkdownViewer from 'etna-js/components/markdown_viewer';

export const MarkdownName = (props) => {
  let {item} = props;
  let {attribute_name} = item;

  if (attribute_name) return <MagmaName {...props} />;

  return <DefaultName {...props} />;
};

const MarkdownItem = ({item, ...props}) => {
  let {attribute_name, text} = item;

  if (attribute_name) {
    return (
      <div className='item'>
        <AttributeViewer
          {...props}
          component={MarkdownAttribute}
          attribute_name={attribute_name}
        />
      </div>
    );
  }

  return (
    <div className='item'>
      <MarkdownViewer text={text} />
    </div>
  );
};

export default MarkdownItem;
