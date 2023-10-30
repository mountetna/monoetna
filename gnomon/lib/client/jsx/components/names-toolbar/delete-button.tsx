import React from "react";
import { useSelector} from 'react-redux';
import DeleteOutlineOutlinedIcon from "@material-ui/icons/DeleteOutlineOutlined";

import { useDispatch } from "../../utils/redux";
import ToolbarButtonWithPopper from "./toolbar-button-with-popper";
import { deleteSelectedGroupsWithNames } from "../../actions/names";
import { selectGlobalState } from "../../selectors/global";
import { selectSelectedCreateNameGroupIds } from "../../selectors/names";



const DeleteButton = ({ small }: { small: boolean }) => {
    const dispatch = useDispatch()

    const globalState = useSelector(selectGlobalState)
    const selectedCreateNameGroupLocalIds = useSelector(selectSelectedCreateNameGroupIds)

    const handleClickDelete = () => {
        dispatch(deleteSelectedGroupsWithNames(globalState))
    }

    return (
        <ToolbarButtonWithPopper
            text="Delete"
            iconComponent={<DeleteOutlineOutlinedIcon />}
            variant={small ? "compact" : "full"}
            color="primary"
            onClickOrPopperChange={handleClickDelete}
            disabled={selectedCreateNameGroupLocalIds.size == 0}
        />
    )
};

export default DeleteButton;