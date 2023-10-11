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
import { selectCreateNameGroupsByLocalId, selectCreateNameLocalIdsByGroupId, selectCreateNameTokenValueLocalIdsByCreateNameLocalId, selectCreateNameTokenValuesByLocalId, selectCreateNamesByLocalId, selectSelectedCreateNameGroupIds } from "../../selectors/names"
import { selectCommonRulesFromSelection } from "../../selectors/global";
import { duplicateCreateNameGroups, iterateOnCreateNameGroupsByRule } from "../../actions/names";



const CopyCreateNameGroupRadio = ({ radioValue, label, quantityValue, onChangeQuantity }: {
    radioValue: string,
    label: string,
    quantityValue?: number,
    onChangeQuantity: (value?: number) => void
}) => {

    const handleChangeQuantity = (event: React.ChangeEvent) => {
        const eventValue = event.target.value
        const quantity = eventValue == "" ? undefined : Number(eventValue)
        onChangeQuantity(quantity)
    }

    return (
        <FormGroup className={"dsadsa"}>
            <FormControl>
                <Radio value={radioValue} />
                <FormLabel>{label}</FormLabel>
            </FormControl>
            <FormControl>
                <InputBase
                    value={quantityValue != undefined ? String(quantityValue) : ""}
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
            <FormControl>
                <Radio value={radioValue} />
                <FormLabel>{label}</FormLabel>
            </FormControl>
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


const ReplaceValuesRadio = () => {

}


const addFromSelectionButtonUseStyles = makeStyles((theme) => ({
    addFromSelectionContainer: {
        display: "inline-block",
    },
}));


const CopyAndReplaceButton = () => {
    const classes = addFromSelectionButtonUseStyles()
    const [open, setOpen] = useState<boolean>(false);
    const [rule, setRule] = useState<Rule>()
    const [radioValue, setRadioValue] = useState<string>("copy")
    const [boundaries, setBoundaries] = useState<IterationBoundaries>({ min: undefined, max: undefined })
    const [quantity, setQuantity] = useState<number | undefined>(undefined)
    const dispatch = useDispatch()
    const anchorEl = useRef(null)

    const selectedCreateNameGroupLocalIds: string[] = useSelector(selectSelectedCreateNameGroupIds)
    const createNameGroupsByLocalId: Record<string, CreateNameGroup> = useSelector(selectCreateNameGroupsByLocalId)

    const selectedCreateNameGroups: CreateNameGroup[] = [...selectedCreateNameGroupLocalIds].map(cngLocalId => createNameGroupsByLocalId[cngLocalId])

    const createNameLocalIdsByCreateNameGroupLocalId: Record<string, string[]> = useSelector(selectCreateNameLocalIdsByGroupId)
    const createNameTokenValueLocalIdsByCreateNameLocalId: Record<string, string[]> = useSelector(selectCreateNameTokenValueLocalIdsByCreateNameLocalId)
    const createNameTokenValuesByLocalId: Record<string, CreateNameTokenValue> = useSelector(selectCreateNameTokenValuesByLocalId)
    const createNamesByLocalId: Record<string, CreateName> = useSelector(selectCreateNamesByLocalId)

    const iterableRules: Rule[] = useSelector(selectCommonRulesFromSelection).filter(rule => rule.hasCounter)

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
                dispatch(duplicateCreateNameGroups(
                    selectedCreateNameGroups,
                    createNameLocalIdsByCreateNameGroupLocalId,
                    createNamesByLocalId,
                    createNameTokenValueLocalIdsByCreateNameLocalId,
                    createNameTokenValuesByLocalId,
                    quantity,
                    {},
                ))
                break
            case "iterate":
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
                break
            default:
                console.error(`Unsupported radio value: ${radioValue}`)
        }

        handleClose()
    }

    return (
        <div className={classes.addFromSelectionContainer}>
            <Button
                startIcon={<LibraryAddOutlinedIcon />}
                onClick={handleToggle}
                ref={anchorEl}
                color="primary"
                aria-label="Copy and Replace"
                aria-haspopup="true"
                aria-controls={open ? "copy-and-replace-dialogue" : undefined}
                disableRipple
                disableFocusRipple
                disableElevation
                disabled={selectedCreateNameGroupLocalIds.length == 0}
            >
                Copy and Replace
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
                                <FormControl component="fieldset" id="copy-and-replace-dialogue">
                                    <RadioGroup
                                        aria-label="Copy and Replace Options"
                                        name="addFromSelectionOptions"
                                        value={radioValue}
                                        onChange={handleChangeRadioValue}
                                    >
                                        <CopyCreateNameGroupRadio
                                            radioValue="copy"
                                            label="How many?"
                                            quantityValue={quantity}
                                            onChangeQuantity={setQuantity}
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

export default CopyAndReplaceButton;