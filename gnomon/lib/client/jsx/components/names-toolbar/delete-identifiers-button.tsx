import React, { useState, useRef } from 'react';
import { useAppSelector as useSelector } from '../../hooks';
import DeleteOutlineOutlinedIcon from '@material-ui/icons/DeleteOutlineOutlined';

import { useDispatch } from '../../utils/redux';
import ToolbarButtonWithPopper from './toolbar-button-with-popper';
import ConfirmationPopper from '../confirmation-popper';
import MessageModal from '../message-modal';
import ConfirmationCodeModal from '../confirmation-code-modal';
import { json_get, json_post } from 'etna-js/utils/fetch';
import { magmaPath } from 'etna-js/api/magma_api';



const DeleteIdentifiersButton = ({ small, project_name, data, buttonText, className, refresh }: { small: boolean, project_name: string, className?: string, data: Array<any>, buttonText: string, refresh: Function }) => {
    const [confirmation, setConfirmation] = useState('');
    const [confirmationError, setConfirmationError] = useState('');
    const [ message, setMessage ] = useState('');
    const [ messageTitle, setMessageTitle ] = useState('');

    const handleDelete = () => {
      json_post(magmaPath(`gnomon/${project_name}/delete`),
      { identifiers: data.map(d => d.name), confirmation }).then(
        ({success}) => {
          setMessage(success);
          setMessageTitle('Deleted Identifiers');
          setConfirmation('');
          refresh();
        }
      ).catch(
        (error) => error.then(
          ({error}:{ error: string}) => {
            setConfirmation('');
            if (error.includes('Confirm')) setConfirmationError(error);
            else {
              setMessage(error.length > 100 ? `${error.substring(0,100)}...` : error);
              setMessageTitle('Deletion Failed');
            }
          }
        )
      )
    };

    return (
      <React.Fragment>
        <ToolbarButtonWithPopper
          text={ buttonText }
          iconComponent={<DeleteOutlineOutlinedIcon />}
          variant={small ? 'compact' : 'full'}
          color="#FF0000"
          onClickOrPopperChange={handleDelete}
          disabled={data.length == 0}
          className={className}
        />
        <ConfirmationCodeModal
          open={!!confirmationError}
          confirmationError={ confirmationError }
          confirmation={ confirmation }
          onChange={ setConfirmation }
          onClose={ () => {
            setConfirmationError('');
            setConfirmation('');
          } }
          onConfirm={ () => {
            setConfirmationError('');
            handleDelete();
          }}
        />
        <MessageModal title={messageTitle} open={!!message} onClose={ () => setMessage('') } message={message}/>
      </React.Fragment>
    );
};

export default DeleteIdentifiersButton;
