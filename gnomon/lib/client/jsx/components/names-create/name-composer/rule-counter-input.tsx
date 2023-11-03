import React from "react";
import { makeStyles } from '@material-ui/core/styles';
import { useSelector } from 'react-redux'
import ButtonBase from "@material-ui/core/ButtonBase";
import FormHelperText from '@material-ui/core/FormHelperText';
import FormControl from '@material-ui/core/FormControl';
import _ from "lodash";

import { selectCompleteCreateNameParentLocalIdsByRenderedValues, selectRenderedCompleteCreateNamesByLocalId } from "../../../selectors/names";
import { fetchNextCounterValueFromMagma } from "../../../utils/names";
import { selectRuleParentLocalIdsByRuleName } from "../../../selectors/rules";
import { UNSET_VALUE } from "../../../models";
import AutosizeTextInput from "../../../utils/AutosizeInput";



const useStyles = makeStyles((theme) => {
    const fontSize = "16px"
    const lineHeight = "1em"

    return {
        ruleCounterField: {
            display: "inline-block",
            lineHeight: lineHeight,
            "&:not(:last-child)": {
                marginRight: "0.35em",
            },
        },
        form: {
            alignItems: "flex-end",
        },
        autoIncrementButton: {
            opacity: "0.25",
            transition: "opacity 0.2s ease-in",
            position: "absolute",
            bottom: "1.5em",
            fontSize: fontSize,
            lineHeight: lineHeight,
            "&:hover, &:active": {
                opacity: "1"
            },
        },
        ruleCounterInput: {
            fontWeight: "bold",
            fontSize: fontSize,
            "& input": {
                "-moz-appearance": "textfield",
                padding: "0",
                height: "unset",
                lineHeight: lineHeight,
            },
            "& input::-webkit-outer-spin-button, & input::-webkit-inner-spin-button": {
                display: "none",
                margin: "0",
            },
        },
        unset: {
            color: "red",
            "& button": {
                opacity: "1",
            },
            "& input::placeholder": {
                color: "red",
            },
        },
    }
});


const RuleCounterField = ({
    value,
    renderedTokensPrefix,
    parentCompleteCreateNameLocalId,
    projectName,
    ruleName,
    handleSetCounterValue
}: {
    value: number | undefined,
    renderedTokensPrefix: string | undefined,
    // TODO: handle multiple parents
    parentCompleteCreateNameLocalId: string | undefined,
    projectName: string,
    ruleName: string,
    handleSetCounterValue: (value?: number) => void,
}) => {

    const classes = useStyles()

    const completeCreateNameParentLocalIdsByRenderedValues = useSelector(selectCompleteCreateNameParentLocalIdsByRenderedValues)
    const needsParentCompleteCreateName = useSelector(selectRuleParentLocalIdsByRuleName)[ruleName] != undefined
    let fullRenderedTokensPrefix = useSelector(selectRenderedCompleteCreateNamesByLocalId)[parentCompleteCreateNameLocalId]

    if (
        renderedTokensPrefix != undefined
        && parentCompleteCreateNameLocalId != undefined
    ) {
        fullRenderedTokensPrefix += renderedTokensPrefix

    } else if (!needsParentCompleteCreateName) {
        fullRenderedTokensPrefix = renderedTokensPrefix
    }

    const hasValue = value != undefined

    const handleClickAutoIncrement = () => {
        if (!(
            renderedTokensPrefix != undefined
            && (
                parentCompleteCreateNameLocalId != undefined
                || !needsParentCompleteCreateName
            )
            && fullRenderedTokensPrefix != undefined
        )) { return }

        const ccnpLocalId: string = parentCompleteCreateNameLocalId != undefined ? parentCompleteCreateNameLocalId : UNSET_VALUE
        const hierarchyValues: number[] = []

        if (
            ccnpLocalId in completeCreateNameParentLocalIdsByRenderedValues
            && renderedTokensPrefix in completeCreateNameParentLocalIdsByRenderedValues[ccnpLocalId]
        ) {
            hierarchyValues.push(
                ...Object.keys(completeCreateNameParentLocalIdsByRenderedValues[ccnpLocalId][renderedTokensPrefix])
                    .map(el => {
                        const asNum = Number.parseInt(el)
                        return Number.isNaN(asNum) ? -1 : asNum
                    })
            )
        }

        const localMaxValue = hierarchyValues.length ? Math.max(...hierarchyValues) : -1

        fetchNextCounterValueFromMagma(projectName, ruleName, fullRenderedTokensPrefix)
            .then(remoteNextValue => {
                handleSetCounterValue(Math.max(localMaxValue + 1, remoteNextValue))
            })
            .catch(err => {
                console.error(`Error auto-incrementing counterValue: ${err}`)
            })
    }

    const handleChangeInput = (event: React.ChangeEvent<HTMLInputElement>) => {
        const eventValue = event.target.value
        const counterValue = eventValue == "" ? undefined : Number(eventValue)
        handleSetCounterValue(counterValue)
    }

    return (
        <span className={classes.ruleCounterField + " " + (!hasValue ? `${classes.unset}` : "")}>
            <FormControl className={classes.form}>
                <ButtonBase
                    onClick={handleClickAutoIncrement}
                    aria-label={`Set Counter Value for "${ruleName}"`}
                    aria-haspopup="false"
                    aria-controls="rule-counter-input"
                    disableRipple
                    disableTouchRipple
                    className={classes.autoIncrementButton}
                >
                    +1
                </ButtonBase>
                <AutosizeTextInput
                    value={hasValue ? String(value) : ""}
                    onChange={handleChangeInput}
                    type="number"
                    inputMode="numeric"
                    inputProps={{ min: 0 }}
                    id="rule-counter-input"
                    placeholder="n"
                    className={classes.ruleCounterInput}
                    minWidth="0.65em"
                />
            </FormControl>
        </span>
    )
}


export default RuleCounterField