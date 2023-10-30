import React, { useState, useEffect } from "react";
import { useSelector, batch } from 'react-redux'
import Button from "@material-ui/core/Button";
import LibraryAddOutlinedIcon from "@material-ui/icons/LibraryAddOutlined";
import Radio from '@material-ui/core/Radio';
import RadioGroup from '@material-ui/core/RadioGroup';
import FormControlLabel from '@material-ui/core/FormControlLabel';
import FormGroup from '@material-ui/core/FormGroup';
import FormLabel from '@material-ui/core/FormLabel';
import FormControl from '@material-ui/core/FormControl';
import InputBase from '@material-ui/core/InputBase';
import Grid from '@material-ui/core/Grid';
import _ from "lodash"

import { CreateNameGroup, Rule } from "../../models";
import { RuleSelect } from "../names-create/name-composer/select";
import CreateNameGroupComposer from "../names-create/name-composer/name-composer";
import { inputEventValueToNumber } from "../../utils/input";
import { selectCreateNameGroupsByLocalId, selectSelectedCreateNameGroupIds } from "../../selectors/names"
import { selectCommonRulesFromSelection, selectGlobalState, selectReplaceRuleFromSelection } from "../../selectors/global";
import { addCreateNameGroupsToReplaceCriteria, createNamesWithGroupForRule, deleteGroupsWithNames, duplicateCreateNameGroups, iterateOnCreateNameGroupsByRule, addOrReplaceCreateNameTokenValuesAndRuleCounterValuesFromReplaceCriteria } from "../../actions/names";
import { useDispatch } from "../../utils/redux";
import ToolbarButtonWithPopper from "./toolbar-button-with-popper";



const CopyCreateNameGroupRadio = ({ radioValue, label, quantityValue, onChangeQuantity }: {
    radioValue: string,
    label: string,
    quantityValue: number | undefined,
    onChangeQuantity: (value?: number) => void
}) => {

    const handleChangeQuantity = (event: React.ChangeEvent<HTMLInputElement>) => {
        const eventValue = event.target.value
        const quantity = eventValue == "" ? undefined : Number(eventValue)
        onChangeQuantity(quantity)
    }

    return (
        <FormGroup className={"dsadsa"}>
            <FormControl>
                <Radio value={radioValue} disableRipple />
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
    start: number | undefined
    end: number | undefined
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

    const handleChangeBoundaries = (start: number | undefined, end: number | undefined) => {
        onChangeBoundaries({ start, end })
    }

    const getBoundaryValue = (boundary: number | undefined): string => {
        return boundary != undefined ? String(boundary) : ""
    }

    return (
        <FormGroup>
            <FormControl>
                <Radio value={radioValue} disableRipple />
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
                    onChange={(event: React.ChangeEvent<HTMLInputElement>) => handleChangeBoundaries(inputEventValueToNumber(event), boundaries.end)}
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
                    onChange={(event: React.ChangeEvent<HTMLInputElement>) => handleChangeBoundaries(boundaries.start, inputEventValueToNumber(event))}
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


const ReplaceValuesForRuleRadio = ({ radioValue, label, createNameGroup }: {
    radioValue: string,
    label: string,
    createNameGroup?: CreateNameGroup,
}) => {
    return (
        <FormGroup>
            <FormControl>
                <Radio value={radioValue} disableRipple />
                <FormLabel>{label}</FormLabel>
            </FormControl>
            {
                createNameGroup ?
                    <CreateNameGroupComposer
                        createNameGroup={createNameGroup}
                    /> : "Select a Name first"
            }
        </FormGroup>
    )
}


const CopyAndReplaceButton = ({ small }: { small: boolean }) => {
    const [open, setOpen] = useState<boolean>(false);
    const [rule, setRule] = useState<Rule>()
    const [radioValue, setRadioValue] = useState<string>("copy")
    const [boundaries, setBoundaries] = useState<IterationBoundaries>({ start: undefined, end: undefined })
    const [quantity, setQuantity] = useState<number | undefined>(undefined)
    const [replaceCreateNameGroup, setReplaceCreateNameGroup] = useState<CreateNameGroup | undefined>()

    const dispatch = useDispatch()

    const globalState = useSelector(selectGlobalState)
    const selectedCreateNameGroupLocalIds = useSelector(selectSelectedCreateNameGroupIds)
    const createNameGroupsByLocalId = useSelector(selectCreateNameGroupsByLocalId)

    const selectedCreateNameGroups = [...selectedCreateNameGroupLocalIds].map(cngLocalId => createNameGroupsByLocalId[cngLocalId])

    const iterableRules = useSelector(selectCommonRulesFromSelection).filter(rule => rule.hasCounter)

    const replaceRule = useSelector(selectReplaceRuleFromSelection)

    // manage CreateNameGroup for <ReplaceValuesForRuleRadio />
    useEffect(() => {
        if (!replaceRule) {
            if (replaceCreateNameGroup) {
                batch(() => {
                    dispatch(deleteGroupsWithNames([replaceCreateNameGroup.localId], globalState))
                    setReplaceCreateNameGroup(undefined)
                })
                return
            }
            setReplaceCreateNameGroup(undefined)
            return
        }

        const actionPayload = createNamesWithGroupForRule(replaceRule.name, globalState, false)
        const newCng = actionPayload.createNameGroups[0]

        if (replaceCreateNameGroup) {
            batch(() => {
                dispatch(deleteGroupsWithNames([replaceCreateNameGroup.localId], globalState))
                dispatch(actionPayload)
                dispatch(addCreateNameGroupsToReplaceCriteria([newCng.localId]))
                setReplaceCreateNameGroup(newCng)
            })
        } else {
            batch(() => {
                dispatch(actionPayload)
                dispatch(addCreateNameGroupsToReplaceCriteria([newCng.localId]))
                setReplaceCreateNameGroup(newCng)
            })
        }
    }, [replaceRule?.name])

    const handleToggle = () => {
        setOpen(prev => !prev);
    };

    const handleClose = () => {
        setOpen(false);
    };

    const handleChangeRadioValue = (event: React.ChangeEvent<HTMLInputElement>) => {
        setRadioValue(event.target.value)
    }

    const handleClickAdd = () => {
        switch (radioValue) {
            case "copy":
                dispatch(duplicateCreateNameGroups(
                    selectedCreateNameGroups,
                    globalState,
                    true,
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
                    globalState,
                ))
                break
            case "replace":
                dispatch(addOrReplaceCreateNameTokenValuesAndRuleCounterValuesFromReplaceCriteria(globalState))
                break
            default:
                console.error(`Unsupported radio value: ${radioValue}`)
        }
    }

    return (
        <ToolbarButtonWithPopper
            text="Copy and Replace"
            iconComponent={<LibraryAddOutlinedIcon />}
            variant={small ? "compact" : "full"}
            color="primary"
            popperComponent={
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
                        <ReplaceValuesForRuleRadio
                            radioValue="replace"
                            label="Replace Values"
                            createNameGroup={replaceCreateNameGroup}
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
                                {_.capitalize(radioValue)}
                            </Button>
                        </Grid>
                    </Grid>
                </FormControl>
            }
            popperId="copy-and-replace-dialogue"
            onClickOrPopperChange={(open: boolean) => setOpen(open)}
            popperOpen={open}
            disabled={selectedCreateNameGroupLocalIds.size == 0}
        />
    )
};

export default CopyAndReplaceButton;