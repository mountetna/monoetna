import React, { useState, useRef } from "react";
import { useSelector, useDispatch } from 'react-redux'
import { makeStyles } from '@material-ui/core/styles';
import Button from "@material-ui/core/Button";
import Popper from "@material-ui/core/Popper";
import Grow from "@material-ui/core/Grow";
import ClickAwayListener from "@material-ui/core/ClickAwayListener";
import Paper from "@material-ui/core/Paper";
import LibraryAddOutlinedIcon from "@material-ui/icons/LibraryAddOutlined";
import Radio from '@material-ui/core/Radio';
import RadioGroup from '@material-ui/core/RadioGroup';
import FormControlLabel from '@material-ui/core/FormControlLabel';
import FormGroup from '@material-ui/core/FormGroup';
import FormLabel from '@material-ui/core/FormLabel';
import FormControl from '@material-ui/core/FormControl';
import InputBase from '@material-ui/core/InputBase';
import Grid from '@material-ui/core/Grid';
import * as _ from "lodash"

import { CreateName, CreateNameGroup, CreateNameTokenValue, Rule } from "../../models";
import { RuleSelect } from "../names-create/name-composer/select";
import { inputEventValueToNumber } from "../../utils/input";
import { selectCreateNameGroupsByLocalId, selectCreateNameLocalIdsByGroupId, selectCreateNameTokenValueLocalIdsByCreateNameLocalId, selectCreateNameTokenValuesByLocalId, selectCreateNameWithLocalId, selectCreateNamesByLocalId, selectSelectedCreateNameGroupIds } from "../../selectors/names"
import { selectIterableRules } from "../../selectors/global";
import { iterateOnCreateNameGroupsByRule } from "../../actions/names";



const CopyCreateNameGroupRadio = ({ radioValue, label }: { radioValue: string, label: string }) => {
    const [quantity, setQuantity] = useState<number | undefined>(undefined)

    const handleChangeQuantity = (event: React.ChangeEvent) => {
        const eventValue = event.target.value
        const quantityValue = eventValue == "" ? undefined : Number(eventValue)
        setQuantity(quantityValue)
    }

    return (
        <FormGroup className={"dsadsa"}>
            <Radio value={radioValue} />
            <FormControl>
                <FormLabel>{label}</FormLabel>
                <InputBase
                    value={quantity != undefined ? String(quantity) : ""}
                    onChange={handleChangeQuantity}
                    type="number"
                    inputProps={{ min: 1 }}
                    id="copy-radio-quantity"
                    placeholder="n"
                    className={"dsadasdas"}
                />
            </FormControl>
        </FormGroup>
    );
}


interface IterationBoundaries {
    start?: number
    end?: number
}


const IterateRuleRadio = ({ radioValue, label, rules, ruleValue, onChangeRule, boundaries, onChangeBoundaries }: {
    radioValue: string,
    label: string,
    rules: Rule[],
    ruleValue: Rule,
    onChangeRule: (rule: Rule) => void,
    boundaries: IterationBoundaries,
    onChangeBoundaries: (boundaries: IterationBoundaries) => void,
}) => {

    const handleChangeBoundaries = (start?: number, end?: number) => {
        onChangeBoundaries({ start, end })
    }

    const getBoundaryValue = (boundary?: number): string => {
        return boundary != undefined ? String(boundary) : ""
    }

    return (
        <FormGroup>
            <Radio value={radioValue} />
            <FormLabel>{label}</FormLabel>
            <FormControl>
                <RuleSelect
                    values={rules}
                    value={ruleValue}
                    label="Rule"
                    placeholder="rule"
                    onSetRule={onChangeRule}
                />
            </FormControl>
            <FormControl>
                <FormLabel>Start</FormLabel>
                <InputBase
                    value={getBoundaryValue(boundaries.start)}
                    onChange={(event: React.ChangeEvent) => handleChangeBoundaries(inputEventValueToNumber(event), boundaries.end)}
                    type="number"
                    inputProps={{ min: 0 }}
                    id="iterate-radio-boundary-minimum"
                    placeholder="n"
                    className={"dsadasdas"}
                />
            </FormControl>
            <FormControl>
                <FormLabel>End</FormLabel>
                <InputBase
                    value={getBoundaryValue(boundaries.end)}
                    onChange={(event: React.ChangeEvent) => handleChangeBoundaries(boundaries.start, inputEventValueToNumber(event))}
                    type="number"
                    inputProps={{ min: 0 }}
                    id="iterate-radio-boundary-maximum"
                    placeholder="n"
                    className={"dsadasdas"}
                />
            </FormControl>
        </FormGroup>
    )
}


const addFromSelectionButtonUseStyles = makeStyles((theme) => ({
    addFromSelectionContainer: {
        display: "inline-block",
    },
}));


// TODO: split duplicate and iterate into separate buttons
const AddFromSelectionButton = () => {
    const classes = addFromSelectionButtonUseStyles()
    const [open, setOpen] = useState<boolean>(false);
    const [rule, setRule] = useState<Rule>()
    const [radioValue, setRadioValue] = useState<string>("copy")
    const [boundaries, setBoundaries] = useState<IterationBoundaries>({ min: undefined, max: undefined })
    const dispatch = useDispatch()
    const anchorEl = useRef(null)

    const selectedCreateNameGroupLocalIds: string[] = useSelector(selectSelectedCreateNameGroupIds)
    const createNameGroupsByLocalId: Record<string, CreateNameGroup> = useSelector(selectCreateNameGroupsByLocalId)

    const selectedCreateNameGroups: CreateNameGroup[] = selectedCreateNameGroupLocalIds.map(cngLocalId => createNameGroupsByLocalId[cngLocalId])

    const createNameLocalIdsByCreateNameGroupLocalId: Record<string, string[]> = useSelector(selectCreateNameLocalIdsByGroupId)
    const createNameTokenValueLocalIdsByCreateNameLocalId: Record<string, string[]> = useSelector(selectCreateNameTokenValueLocalIdsByCreateNameLocalId)
    const createNameTokenValuesByLocalId: Record<string, CreateNameTokenValue> = useSelector(selectCreateNameTokenValuesByLocalId)
    const createNamesByLocalId: Record<string, CreateName> = useSelector(selectCreateNamesByLocalId)

    const iterableRules: Rule[] = useSelector(selectIterableRules)

    const handleToggle = () => {
        setOpen(prev => !prev);
    };

    const handleClose = () => {
        setOpen(false);
    };

    const handleChangeRadioValue = (event: React.ChangeEvent) => {
        setRadioValue(event.target.value)
    }

    const handleClickAdd = () => {
        switch (radioValue) {
            case "copy":
                console.log(radioValue)
                return
            case "iterate":
                console.log(radioValue)
                dispatch(iterateOnCreateNameGroupsByRule(
                    selectedCreateNameGroups,
                    rule.name,
                    boundaries.start,
                    boundaries.end,
                    createNameLocalIdsByCreateNameGroupLocalId,
                    createNamesByLocalId,
                    createNameTokenValueLocalIdsByCreateNameLocalId,
                    createNameTokenValuesByLocalId,
                ))
                handleClose()
                return
            default:
                console.error(`Unsupported radio value: ${radioValue}`)
                return
        }
    }

    return (
        <div className={classes.addFromSelectionContainer}>
            <Button
                startIcon={<LibraryAddOutlinedIcon />}
                onClick={handleToggle}
                ref={anchorEl}
                color="primary"
                aria-label="Add from Selection"
                aria-haspopup="true"
                aria-controls={open ? "add-from-selection-dialogue" : undefined}
                disableRipple
                disableFocusRipple
                disableElevation
                disabled={selectedCreateNameGroupLocalIds.length == 0}
            >
                Add from Selection
            </Button>
            <Popper
                open={open}
                anchorEl={anchorEl.current}
                placement="bottom"
                role={undefined}
                transition
            >
                {({ TransitionProps }) => (
                    <Grow
                        {...TransitionProps}
                        style={{ transformOrigin: "center top" }}
                    >
                        <Paper variant="outlined">
                            <ClickAwayListener onClickAway={handleClose}>
                                <FormControl component="fieldset">
                                    <RadioGroup
                                        aria-label="Add from Selection Options"
                                        name="addFromSelectionOptions"
                                        value={radioValue}
                                        onChange={handleChangeRadioValue}
                                    >
                                        <CopyCreateNameGroupRadio
                                            radioValue="copy"
                                            label="How many?"
                                        />
                                        <IterateRuleRadio
                                            radioValue="iterate"
                                            label="Iterate on Rule"
                                            rules={iterableRules}
                                            ruleValue={rule}
                                            onChangeRule={setRule}
                                            boundaries={boundaries}
                                            onChangeBoundaries={setBoundaries}
                                        />
                                    </RadioGroup>
                                    <Grid container>
                                        <Grid item xs={6}>
                                            <Button
                                                onClick={handleClose}
                                                color="secondary"
                                                disableElevation
                                            >
                                                Cancel
                                            </Button>
                                        </Grid>
                                        <Grid item xs={6}>
                                            <Button
                                                onClick={handleClickAdd}
                                                color="primary"
                                                disableElevation
                                            >
                                                Add
                                            </Button>
                                        </Grid>
                                    </Grid>
                                </FormControl>
                            </ClickAwayListener>
                        </Paper>
                    </Grow>
                )}
            </Popper>
        </div>
    )
};

export default AddFromSelectionButton;