import React from "react";
import { makeStyles } from '@material-ui/core/styles';
import { useSelector } from 'react-redux'
import ButtonBase from "@material-ui/core/ButtonBase";
import InputBase from "@material-ui/core/InputBase";
import FormHelperText from '@material-ui/core/FormHelperText';
import FormControl from '@material-ui/core/FormControl';
import * as _ from "lodash"

import { CompleteCreateNameParentsByRenderedValues } from "../../../reducers/names";
import { selectCompleteCreateNamesByParentAndValues, selectRenderedCompleteCreateNameWithLocalId } from "../../../selectors/names";
import { fetchNextCounterValueFromMagma } from "../../../utils/names";
import { selectRuleParentLocalIdsByRuleName } from "../../../selectors/rules";
import { UNSET_VALUE } from "../../../models";



const useStyles = makeStyles((theme) => ({
    ruleCounterField: {
        display: "inline-flex",
        flexDirection: "column",
        alignItems: "center",
    },
    ruleCounterInput: (inputWidthEm: number) => ({
        fontWeight: "bold",
        width: `${inputWidthEm}em`,
        minWidth: "1.5em",
        maxWidth: "4em",
    }),
    unsetOrError: {
        color: "red"
    },
    autoIncrementButton: {},
}));


const RuleCounterField = ({
    value,
    renderedTokensPrefix,
    completeCreateNameParentLocalId,
    projectName,
    ruleName,
    handleSetCounterValue
}: {
    value: number | undefined,
    renderedTokensPrefix: string | undefined,
    completeCreateNameParentLocalId: string | undefined,
    projectName: string,
    ruleName: string,
    handleSetCounterValue: (value?: number) => void,
}) => {

    const classes = useStyles(value ? String(value).length : 1)
    const completeCreateNamesByParentAndValues: CompleteCreateNameParentsByRenderedValues = useSelector(selectCompleteCreateNamesByParentAndValues)
    const needsParentCompletedCreateName: boolean = useSelector(selectRuleParentLocalIdsByRuleName)[ruleName] != undefined

    let fullRenderedTokensPrefix: string | undefined = undefined
    if (
        renderedTokensPrefix != undefined
        && completeCreateNameParentLocalId != undefined
    ) {
        fullRenderedTokensPrefix = useSelector(
            state => selectRenderedCompleteCreateNameWithLocalId(state, completeCreateNameParentLocalId)
        )
        fullRenderedTokensPrefix += renderedTokensPrefix

    } else if (!needsParentCompletedCreateName) {
        fullRenderedTokensPrefix = renderedTokensPrefix
    }

    const hasValue = value != undefined
    let counterValueCollision = false  // TODO

    const handleClickAutoIncrement = () => {
        if (!(
            renderedTokensPrefix != undefined
            && (
                completeCreateNameParentLocalId != undefined
                || !needsParentCompletedCreateName
            )
            && fullRenderedTokensPrefix != undefined
        )) { return }

        const ccnpLocalId: string = completeCreateNameParentLocalId != undefined ? completeCreateNameParentLocalId : UNSET_VALUE

        const localMaxValue = ccnpLocalId in completeCreateNamesByParentAndValues
            ? Math.max(
                ...Object.keys(completeCreateNamesByParentAndValues[ccnpLocalId])
                    .map(el => {
                        const asNum = Number.parseInt(el)
                        return Number.isNaN(asNum) ? -1 : asNum
                    })
            )
            : -1

        fetchNextCounterValueFromMagma(projectName, ruleName, fullRenderedTokensPrefix)
            .then(remoteMaxValue => {
                handleSetCounterValue(Math.max(localMaxValue, remoteMaxValue) + 1)
            })
            .catch(err => {
                console.error(`Error auto-incrementing counterValue: ${err}`)
            })
    }

    const handleChangeInput = (event: React.ChangeEvent) => {
        const eventValue = event.target.value
        const counterValue = eventValue == "" ? undefined : Number(eventValue)
        handleSetCounterValue(counterValue)
    }

    return (
        <span className={classes.ruleCounterField + " " + ((!hasValue || counterValueCollision) ? `${classes.unsetOrError}` : "")}>
            <FormControl error={counterValueCollision}>
                <ButtonBase
                    onClick={handleClickAutoIncrement}
                    aria-label={`Set Counter Value for "${ruleName}"`}
                    aria-haspopup="false"
                    aria-controls="rule-counter-input"
                    disableRipple
                    disableTouchRipple
                    className={classes.autoIncrementButton}
                    disabled={hasValue && !counterValueCollision}
                >
                    +1
                </ButtonBase>
                <InputBase
                    value={hasValue ? String(value) : ""}
                    onChange={handleChangeInput}
                    type="number"
                    inputProps={{ min: 0 }}
                    id="rule-counter-input"
                    placeholder="n"
                    className={classes.ruleCounterInput}
                />
                {/* TODO: make this a tooltip? */}
                <FormHelperText>
                    {counterValueCollision && "Value already exists"}
                </FormHelperText>
            </FormControl>
        </span>
    )
}


export default RuleCounterField