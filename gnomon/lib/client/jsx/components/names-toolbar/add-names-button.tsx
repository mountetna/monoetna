import React, { useState } from "react";
import { useSelector, batch } from 'react-redux';
import MenuList from "@material-ui/core/MenuList";
import MenuItem from "@material-ui/core/MenuItem";
import AddCircleOutlineIcon from "@material-ui/icons/AddCircleOutline";

import { selectRulesByName } from "../../selectors/rules";
import { createNamesWithGroupForRule } from "../../actions/names";
import { useDispatch } from "../../utils/redux";
import { selectGlobalState } from "../../selectors/global";
import ToolbarButtonWithPopper from "./toolbar-button-with-popper";



const AddNamesButton = ({ small }: { small: boolean }) => {
    const [open, setOpen] = useState<boolean>(false);
    const dispatch = useDispatch()

    const rules = useSelector(selectRulesByName)
    const globalState = useSelector(selectGlobalState)

    const handleClose = () => {
        setOpen(false);
    };

    const handleClickRule = (ruleName: string) => {
        batch(() => {
            dispatch(createNamesWithGroupForRule(ruleName, globalState, true))
            handleClose();
        })
    };

    return (
        <ToolbarButtonWithPopper
            text="Add Names"
            iconComponent={<AddCircleOutlineIcon />}
            variant={small ? "compact" : "full"}
            color="primary"
            popperComponent={
                <MenuList autoFocusItem={open} id="add-names-rules-menu">
                    {
                        Object.keys(rules).map((ruleName) =>
                            <MenuItem
                                onClick={() => handleClickRule(ruleName)}
                                key={ruleName}
                                disableRipple
                            >
                                {ruleName}
                            </MenuItem>
                        )
                    }
                </MenuList>
            }
            popperId="add-names-rules-menu"
            onClickOrPopperChange={(open: boolean) => setOpen(open)}
            popperOpen={open}
        />
    )
};

export default AddNamesButton;