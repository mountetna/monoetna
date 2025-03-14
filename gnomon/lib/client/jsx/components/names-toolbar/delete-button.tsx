import React, { useState, useRef } from 'react';
import { batch } from 'react-redux';
import { useAppSelector as useSelector } from '../../hooks';
import DeleteOutlineOutlinedIcon from '@material-ui/icons/DeleteOutlineOutlined';

import { useDispatch } from '../../utils/redux';
import ToolbarButtonWithPopper from './toolbar-button-with-popper';
import { deleteSelectedGroupsWithNames } from '../../actions/names';
import { selectGlobalState } from '../../selectors/global';
import { selectSelectedCreateNameGroupIds } from '../../selectors/names';
import ConfirmationPopper from '../confirmation-popper';



const DeleteButton = ({ small, className }: { small: boolean, className?: string }) => {
    const dispatch = useDispatch();

    const [confirmationOpen, setConfirmationOpen] = useState<boolean>(false);

    const buttonRef = useRef(null);

    const globalState = useSelector(selectGlobalState);
    const selecedCount = useSelector(selectSelectedCreateNameGroupIds).size;

    const handleConfirmDelete = (confirmed: boolean) => {
        if (!confirmed) {
            setConfirmationOpen(false);
            return;
        }

        batch(() => {
            dispatch(deleteSelectedGroupsWithNames(globalState));
            setConfirmationOpen(false);
        });
    };

    return (
        <React.Fragment>
            <ToolbarButtonWithPopper
                text="Delete Selection"
                iconComponent={<DeleteOutlineOutlinedIcon />}
                variant={small ? 'compact' : 'full'}
                color="#FF0000"
                onClickOrPopperChange={() => setConfirmationOpen(open => !open)}
                disabled={selecedCount == 0}
                buttonRef={buttonRef}
                className={className}
            />
            <ConfirmationPopper
                open={confirmationOpen}
                onConfirm={handleConfirmDelete}
                onClose={() => setConfirmationOpen(false)}
                anchorRef={buttonRef}
            />
        </React.Fragment>
    );
};

export default DeleteButton;