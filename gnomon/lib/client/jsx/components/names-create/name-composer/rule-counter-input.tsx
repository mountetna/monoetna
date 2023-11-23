import React from 'react';
import { makeStyles } from '@material-ui/core/styles';
import { useSelector } from 'react-redux';
import ButtonBase from '@material-ui/core/ButtonBase';
import FormControl from '@material-ui/core/FormControl';
import _ from 'lodash';

import { selectCompleteCreateNameParentLocalIdsByRenderedValues, selectRenderedCompleteCreateNamesByLocalId } from '../../../selectors/names';
import { fetchNextCounterValueFromMagma } from '../../../utils/names';
import { selectRuleParentLocalIdsByRuleName } from '../../../selectors/rules';
import { UNSET_VALUE } from '../../../models';
import AutosizeTextInput from '../../autosize-text-input';



const useStyles = makeStyles((theme) => {
    const fontSize = '16px';
    const lineHeight = '1em';

    return {
        ruleCounterField: {
            display: 'inline-block',
            lineHeight: lineHeight,
            '&:not(:last-child)': {
                marginRight: '0.35em',
            },
        },
        form: {
            alignItems: 'flex-end',
        },
        autoIncrementButton: {
            opacity: '0.25',
            transition: 'opacity 0.2s ease-in',
            position: 'absolute',
            bottom: '1.5em',
            fontSize: fontSize,
            lineHeight: lineHeight,
            '&:hover, &:active': {
                opacity: '1'
            },
        },
        ruleCounterInput: {
            fontWeight: 'bold',
            fontSize: fontSize,
            background: 'none',
            '&.highlighted': {
                background: 'yellow',
            },
            '& input': {
                '-moz-appearance': 'textfield',
                padding: '0',
                height: 'unset',
                lineHeight: lineHeight,
            },
            '& input::-webkit-outer-spin-button, & input::-webkit-inner-spin-button': {
                display: 'none',
                margin: '0',
            },
        },
        unset: {
            color: 'red',
            '& button': {
                opacity: '1',
            },
            '& input::placeholder': {
                opacity: '1',
                color: 'red',
            },
        },
    };
});


const RuleCounterField = ({
    value,
    renderedTokensPrefix,
    parentCompleteCreateNameLocalId,
    projectName,
    ruleName,
    includeRuleCounterIncrementer = true,
    highlight,
    handleSetCounterValue
}: {
    value: number | undefined,
    renderedTokensPrefix: string | undefined,
    // TODO: handle multiple parents
    parentCompleteCreateNameLocalId: string | undefined,
    projectName: string,
    ruleName: string,
    includeRuleCounterIncrementer?: boolean,
    highlight: boolean,
    handleSetCounterValue: (value?: number) => void,
}) => {

    const classes = useStyles({ highlight });

    const completeCreateNameParentLocalIdsByRenderedValues = useSelector(selectCompleteCreateNameParentLocalIdsByRenderedValues);
    const needsParentCompleteCreateName = useSelector(selectRuleParentLocalIdsByRuleName)[ruleName] != undefined;
    const renderedCompleteCreateNamesByLocalId = useSelector(selectRenderedCompleteCreateNamesByLocalId);

    let fullRenderedTokensPrefix: string | undefined;

    if (parentCompleteCreateNameLocalId) {
        fullRenderedTokensPrefix = renderedCompleteCreateNamesByLocalId[parentCompleteCreateNameLocalId];
    }

    if (
        renderedTokensPrefix != undefined
        && parentCompleteCreateNameLocalId != undefined
    ) {
        fullRenderedTokensPrefix += renderedTokensPrefix;

    } else if (!needsParentCompleteCreateName) {
        fullRenderedTokensPrefix = renderedTokensPrefix;
    }

    const hasValue = value != undefined;

    const handleClickAutoIncrement = async () => {
        if (!(
            renderedTokensPrefix != undefined
            && (
                parentCompleteCreateNameLocalId != undefined
                || !needsParentCompleteCreateName
            )
            && fullRenderedTokensPrefix != undefined
        )) { return; }

        const ccnpLocalId: string = parentCompleteCreateNameLocalId != undefined ? parentCompleteCreateNameLocalId : UNSET_VALUE;
        const hierarchyValues: number[] = [];

        if (
            ccnpLocalId in completeCreateNameParentLocalIdsByRenderedValues
            && renderedTokensPrefix in completeCreateNameParentLocalIdsByRenderedValues[ccnpLocalId]
        ) {
            hierarchyValues.push(
                ...Object.keys(completeCreateNameParentLocalIdsByRenderedValues[ccnpLocalId][renderedTokensPrefix])
                    .map(el => {
                        const asNum = Number.parseInt(el);
                        return Number.isNaN(asNum) ? -1 : asNum;
                    })
            );
        }

        const localMaxValue = hierarchyValues.length ? Math.max(...hierarchyValues) : -1;

        try {
            const remoteNextValue = await fetchNextCounterValueFromMagma(
                projectName,
                ruleName,
                fullRenderedTokensPrefix,
            );
            handleSetCounterValue(Math.max(localMaxValue + 1, remoteNextValue));
        } catch (err) {
            console.error(`Error auto-incrementing counterValue: ${err}`);
        }
    };

    const handleChangeInput = (event: React.ChangeEvent<HTMLInputElement>) => {
        const eventValue = event.target.value;
        const counterValue = eventValue == '' ? undefined : Number(eventValue);
        handleSetCounterValue(counterValue);
    };

    return (
        <span className={classes.ruleCounterField + ' ' + (!hasValue ? `${classes.unset}` : '')}>
            <FormControl className={classes.form}>
                {includeRuleCounterIncrementer && <ButtonBase
                    onClick={handleClickAutoIncrement}
                    aria-label={`Set Counter Value for ${ruleName}`}
                    aria-haspopup="false"
                    disableRipple
                    disableTouchRipple
                    className={classes.autoIncrementButton}
                >
                    +1
                </ButtonBase>}
                <AutosizeTextInput
                    value={hasValue ? String(value) : ''}
                    onChange={handleChangeInput}
                    type="number"
                    inputMode="numeric"
                    inputProps={{ min: 0, 'aria-label': `${ruleName}-counter-value` }}
                    placeholder="n"
                    className={classes.ruleCounterInput + (highlight ? ' highlighted' : '')}
                    minWidth="0.65em"
                />
            </FormControl>
        </span>
    );
};


export default RuleCounterField;