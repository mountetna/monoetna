import React from 'react';
import { batch } from 'react-redux';
import { makeStyles } from '@material-ui/core/styles';
import ButtonBase from '@material-ui/core/ButtonBase';
import FormControl from '@material-ui/core/FormControl';
import _ from 'lodash';

import { selectCompleteCreateNameParentLocalIdsByRenderedValues, selectMagmaIncrementRequestStatusWithCreateNameGroupLocalId, selectRenderedCompleteCreateNamesByLocalId } from '../../../selectors/names';
import { fetchNextCounterValueFromMagma } from '../../../utils/names';
import { selectRuleParentLocalIdsByRuleName } from '../../../selectors/rules';
import { UNSET_VALUE } from '../../../models';
import AutosizeTextInput from '../../autosize-text-input';
import { useAppSelector as useSelector } from '../../../hooks';
import { useDispatch } from '../../../utils/redux';
import { setMagmaIncrementCounterRequest } from '../../../actions/names';



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
                padding: '0',
                height: 'unset',
                lineHeight: lineHeight,
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
    createNameGroupLocalId,
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
    createNameGroupLocalId: string,
    projectName: string,
    ruleName: string,
    includeRuleCounterIncrementer?: boolean,
    highlight: boolean,
    handleSetCounterValue: (value?: number) => void,
}) => {

    const classes = useStyles({ highlight });
    const dispatch = useDispatch();

    const completeCreateNameParentLocalIdsByRenderedValues = useSelector(selectCompleteCreateNameParentLocalIdsByRenderedValues);
    const needsParentCompleteCreateName = useSelector(selectRuleParentLocalIdsByRuleName)[ruleName] != undefined;
    const renderedCompleteCreateNamesByLocalId = useSelector(selectRenderedCompleteCreateNamesByLocalId);

    const incrementRequestStatus = useSelector(state => selectMagmaIncrementRequestStatusWithCreateNameGroupLocalId(state, createNameGroupLocalId));

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

        dispatch(setMagmaIncrementCounterRequest(createNameGroupLocalId, 'inProgress'));

        try {
            const remoteNextValue = await fetchNextCounterValueFromMagma(
                projectName,
                ruleName,
                fullRenderedTokensPrefix,
            );

            batch(() => {
                handleSetCounterValue(Math.max(localMaxValue + 1, remoteNextValue));
                dispatch(setMagmaIncrementCounterRequest(createNameGroupLocalId, 'idle'));
            });
        } catch (err) {
            console.error(`Error auto-incrementing counterValue: ${err}`);
            dispatch(setMagmaIncrementCounterRequest(createNameGroupLocalId, 'error'));
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
                    disabled={incrementRequestStatus && incrementRequestStatus == 'inProgress'}
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