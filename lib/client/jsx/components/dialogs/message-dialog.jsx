import * as React from 'react';
import Icon from '../icon';

const ICONS = {
  notice: 'bell',
  warning: 'exclamation-triangle',
  error: 'skull-crossbones'
};

const MessageDialog = ({title,message,message_type}) =>
   <div className='message-dialog'>
     <div className='title'>{title}</div>
     <div className='message'>
       <Icon icon={ ICONS[message_type] } className={ message_type }/>
       { message }
     </div>
   </div>;

export default MessageDialog;
