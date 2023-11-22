import React, { useState, useEffect, useRef } from 'react';
import { useSelector, batch } from 'react-redux';
import { makeStyles } from '@material-ui/core/styles';
import Button from '@material-ui/core/Button';
import LibraryAddOutlinedIcon from '@material-ui/icons/LibraryAddOutlined';
import Radio from '@material-ui/core/Radio';
import RadioGroup from '@material-ui/core/RadioGroup';
import FormControlLabel from '@material-ui/core/FormControlLabel';
import FormLabel from '@material-ui/core/FormLabel';
import FormControl from '@material-ui/core/FormControl';
import _ from 'lodash';

import { CreateNameGroup, Rule } from '../../models';
import { RuleSelect } from '../names-create/name-composer/select';
import CreateNameGroupComposer from '../names-create/name-composer/name-composer';
import { inputEventValueToNumber } from '../../utils/input';
import { selectCreateNameGroupsByLocalId, selectSelectedCreateNameGroupIds } from '../../selectors/names';
import { selectCommonRulesFromSelection, selectGlobalState, selectReplaceRuleFromSelection } from '../../selectors/global';
import { addCreateNameGroupsToReplaceCriteria, createNamesWithGroupForRule, deleteGroupsWithNames, duplicateCreateNameGroups, iterateOnCreateNameGroupsByRule, addOrReplaceCreateNameTokenValuesAndRuleCounterValuesFromReplaceCriteria } from '../../actions/names';
import { useDispatch } from '../../utils/redux';
import ToolbarButtonWithPopper from './toolbar-button-with-popper';
import AutosizeInput from '../autosize-text-input';



const useStyles = makeStyles((theme) => ({
    formContainer: {
        padding: '1em',
    },
    radioGroup: {
        marginBottom: '2em',
    },
    radioOption: {
        '&:not(:last-child)': {
            marginBottom: '0.5em',
        },
        opacity: '0.5',
        transition: 'opacity 0.2s ease-in',
        '&:hover, &.selected': {
            opacity: '1',
        },
        '& .radio-row': {
            display: 'flex',
            alignItems: 'center',
        },
        '& .form-row': {
            // TODO: make this align perfectly
            marginLeft: '2em',
        },
        '& .form-row > :not(:last-child), & .radio-row > :not(:last-child)': {
            marginRight: '1em',
        },
    },
    iterateBoundaryInputs: {
        flexDirection: 'row',
        alignItems: 'center',
        '& label': {
            marginRight: '1em',
        },
    },
    buttonsContainer: {
        textAlign: 'center',
        '& > button:not(:last-child)': {
            marginRight: '1em',
        },
    },
}));


const CopyCreateNameGroupRadio = ({ radioValue, currentRadioValue, label, quantityValue, onChangeQuantity }: {
    radioValue: string,
    currentRadioValue: string,
    label: string,
    quantityValue: number | undefined,
    onChangeQuantity: (value?: number) => void
}) => {
    const classes = useStyles();

    const handleChangeQuantity = (event: React.ChangeEvent<HTMLInputElement>) => {
        const eventValue = event.target.value;
        const quantity = eventValue == '' ? undefined : Number(eventValue);
        onChangeQuantity(quantity);
    };

    return (
        <div className={`${classes.radioOption} ${radioValue == currentRadioValue ? 'selected' : ''}`}>
            <div className="radio-row">
                <FormControlLabel
                    value={radioValue}
                    label={label}
                    control={<Radio disableRipple />}
                />
                <FormControl>
                    <AutosizeInput
                        value={quantityValue != undefined ? String(quantityValue) : ''}
                        onChange={handleChangeQuantity}
                        type="number"
                        inputProps={{ min: 1, 'aria-label': 'copy-radio-quantity' }}
                        id="copy-radio-quantity"
                        placeholder="n times"
                    />
                </FormControl>
            </div>
        </div>
    );
};


interface IterationBoundaries {
    start: number | undefined
    end: number | undefined
}


const IterateRuleRadio = ({ radioValue, currentRadioValue, label, rules, ruleValue, onChangeRule, boundaries, onChangeBoundaries }: {
    radioValue: string,
    currentRadioValue: string,
    label: string,
    rules: Rule[],
    ruleValue?: Rule,
    onChangeRule: (rule?: Rule) => void,
    boundaries: IterationBoundaries,
    onChangeBoundaries: (boundaries: IterationBoundaries) => void,
}) => {
    const classes = useStyles();

    const minInputRef = useRef();
    const maxInputRef = useRef();

    const handleChangeBoundaries = (start: number | undefined, end: number | undefined) => {
        onChangeBoundaries({ start, end });
    };

    const getBoundaryValue = (boundary: number | undefined): string => {
        return boundary != undefined ? String(boundary) : '';
    };

    const handleFormClick = (event: React.ChangeEvent<HTMLInputElement>) => {
        const container = event.currentTarget;
        container.getElementsByTagName('input')[0].focus();
    };

    return (
        <div className={`${classes.radioOption} ${radioValue == currentRadioValue ? 'selected' : ''}`}>
            <div className="radio-row">
                <FormControlLabel
                    value={radioValue}
                    label={label}
                    control={<Radio disableRipple />}
                />
                <FormControl>
                    <RuleSelect
                        values={rules}
                        value={ruleValue}
                        label="Element"
                        placeholder="element"
                        onSetRule={onChangeRule}
                    />
                </FormControl>
            </div>
            <div className="form-row">
                {/* @ts-ignore */}
                <FormControl className={classes.iterateBoundaryInputs} onClick={handleFormClick}>
                    <FormLabel>Start</FormLabel>
                    <AutosizeInput
                        value={getBoundaryValue(boundaries.start)}
                        onChange={(event: React.ChangeEvent<HTMLInputElement>) => handleChangeBoundaries(inputEventValueToNumber(event), boundaries.end)}
                        type="number"
                        inputProps={{ min: 0, 'aria-label': 'iterate-radio-boundary-minimum' }}
                        id="iterate-radio-boundary-minimum"
                        placeholder="n"
                    />
                </FormControl>
                {/* @ts-ignore */}
                <FormControl className={classes.iterateBoundaryInputs} onClick={handleFormClick}>
                    <FormLabel>End</FormLabel>
                    <AutosizeInput
                        value={getBoundaryValue(boundaries.end)}
                        onChange={(event: React.ChangeEvent<HTMLInputElement>) => handleChangeBoundaries(boundaries.start, inputEventValueToNumber(event))}
                        type="number"
                        inputProps={{ min: 0, 'aria-label': 'iterate-radio-boundary-maximum' }}
                        id="iterate-radio-boundary-maximum"
                        placeholder="n"
                    />
                </FormControl>
            </div>
        </div>
    );
};


const ReplaceValuesForRuleRadio = ({ radioValue, currentRadioValue, label, createNameGroup }: {
    radioValue: string,
    currentRadioValue: string,
    label: string,
    createNameGroup?: CreateNameGroup,
}) => {
    const classes = useStyles();

    return (
        <div className={`${classes.radioOption} ${radioValue == currentRadioValue ? 'selected' : ''}`}>
            <div className="radio-row">
                <FormControlLabel
                    value={radioValue}
                    label={label}
                    control={<Radio disableRipple />}
                />
            </div>
            <div className="form-row">
                {
                    createNameGroup ?
                        <CreateNameGroupComposer
                            createNameGroup={createNameGroup}
                            includeRuleCounterIncrementer={false}
                            checkForDuplicates={false}
                        /> : 'Select a Name first'
                }
            </div>
        </div>
    );
};


const CopyAndReplaceButton = ({ small, className }: { small: boolean, className?: string }) => {
    const classes = useStyles();

    const [open, setOpen] = useState<boolean>(false);
    const [rule, setRule] = useState<Rule>();
    const [radioValue, setRadioValue] = useState<string>('copy');
    const [boundaries, setBoundaries] = useState<IterationBoundaries>({ start: undefined, end: undefined });
    const [quantity, setQuantity] = useState<number | undefined>(undefined);
    const [replaceCreateNameGroup, setReplaceCreateNameGroup] = useState<CreateNameGroup | undefined>();

    const dispatch = useDispatch();

    const globalState = useSelector(selectGlobalState);
    const selectedCreateNameGroupLocalIds = useSelector(selectSelectedCreateNameGroupIds);
    const createNameGroupsByLocalId = useSelector(selectCreateNameGroupsByLocalId);

    const selectedCreateNameGroups = [...selectedCreateNameGroupLocalIds].map(cngLocalId => createNameGroupsByLocalId[cngLocalId]);

    const iterableRules = useSelector(selectCommonRulesFromSelection).filter(rule => rule.hasCounter);

    const replaceRule = useSelector(selectReplaceRuleFromSelection);

    // manage CreateNameGroup for <ReplaceValuesForRuleRadio />
    useEffect(() => {
        if (!replaceRule) {
            if (replaceCreateNameGroup) {
                batch(() => {
                    dispatch(deleteGroupsWithNames([replaceCreateNameGroup.localId], globalState));
                    setReplaceCreateNameGroup(undefined);
                });
                return;
            }
            setReplaceCreateNameGroup(undefined);
            return;
        }

        const actionPayload = createNamesWithGroupForRule(replaceRule.name, globalState, false);
        const newCng = actionPayload.createNameGroups[0];

        if (replaceCreateNameGroup) {
            batch(() => {
                dispatch(deleteGroupsWithNames([replaceCreateNameGroup.localId], globalState));
                dispatch(actionPayload);
                dispatch(addCreateNameGroupsToReplaceCriteria([newCng.localId]));
                setReplaceCreateNameGroup(newCng);
            });
        } else {
            batch(() => {
                dispatch(actionPayload);
                dispatch(addCreateNameGroupsToReplaceCriteria([newCng.localId]));
                setReplaceCreateNameGroup(newCng);
            });
        }
    }, [replaceRule?.name]);

    const handleClose = () => {
        setOpen(false);
    };

    const handleChangeRadioValue = (event: React.ChangeEvent<HTMLInputElement>) => {
        setRadioValue(event.target.value);
    };

    const handleClickAdd = () => {
        switch (radioValue) {
            case 'copy':
                if (quantity == undefined) {
                    break;
                }

                dispatch(duplicateCreateNameGroups(
                    selectedCreateNameGroups,
                    globalState,
                    true,
                    quantity,
                    {},
                ));
                break;
            case 'iterate':
                if (
                    rule == undefined
                    || boundaries.start == undefined
                    || boundaries.end == undefined
                ) {
                    break;
                }

                dispatch(iterateOnCreateNameGroupsByRule(
                    selectedCreateNameGroups,
                    rule.name,
                    boundaries.start,
                    boundaries.end,
                    globalState,
                ));
                break;
            case 'replace':
                dispatch(addOrReplaceCreateNameTokenValuesAndRuleCounterValuesFromReplaceCriteria(globalState));
                break;
            default:
                console.error(`Unsupported radio value: ${radioValue}`);
        }
    };

    return (
        <ToolbarButtonWithPopper
            text="Copy and Replace"
            iconComponent={<LibraryAddOutlinedIcon />}
            variant={small ? 'compact' : 'full'}
            color="primary"
            popperComponent={
                <FormControl
                    component="fieldset"
                    id="copy-and-replace-dialogue"
                    className={classes.formContainer}
                >
                    <RadioGroup
                        aria-label="Copy and Replace Options"
                        name="copy-and-replace-options"
                        value={radioValue}
                        onChange={handleChangeRadioValue}
                        className={classes.radioGroup}
                    >
                        <CopyCreateNameGroupRadio
                            radioValue="copy"
                            currentRadioValue={radioValue}
                            label="Copy"
                            quantityValue={quantity}
                            onChangeQuantity={setQuantity}
                        />
                        <IterateRuleRadio
                            radioValue="iterate"
                            currentRadioValue={radioValue}
                            label="Iterate on Element"
                            rules={iterableRules}
                            ruleValue={rule}
                            onChangeRule={setRule}
                            boundaries={boundaries}
                            onChangeBoundaries={setBoundaries}
                        />
                        <ReplaceValuesForRuleRadio
                            radioValue="replace"
                            currentRadioValue={radioValue}
                            label="Replace Values"
                            createNameGroup={replaceCreateNameGroup}
                        />
                    </RadioGroup>
                    <div className={classes.buttonsContainer}>
                        <Button
                            onClick={handleClose}
                            color="secondary"
                            disableRipple
                            disableElevation
                            aria-label='cancel-copy-or-replace'
                        >
                            Cancel
                        </Button>
                        <Button
                            onClick={handleClickAdd}
                            color="primary"
                            disableRipple
                            disableElevation
                            aria-label='perform-copy-or-replace'
                        >
                            {_.capitalize(radioValue)}
                        </Button>
                    </div>
                </FormControl>
            }
            popperId="copy-and-replace-dialogue"
            onClickOrPopperChange={setOpen}
            popperOpen={open}
            disabled={selectedCreateNameGroupLocalIds.size == 0}
            className={className}
        />
    );
};

export default CopyAndReplaceButton;