import React from "react";
import { makeStyles } from '@material-ui/core/styles';
import { useSelector } from 'react-redux'
import ButtonBase from "@material-ui/core/ButtonBase";
import InputBase from "@material-ui/core/InputBase";
import FormHelperText from '@material-ui/core/FormHelperText';
import FormControl from '@material-ui/core/FormControl';
import * as _ from "lodash"

import { CompleteCreateNameParentsByRenderedValues } from "../../../reducers/names";
import { selectCompleteCreateNameParentLocalIdsByRenderedValues, selectRenderedCompleteCreateNamesByLocalId } from "../../../selectors/names";
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

    const classes = useStyles(value ? String(value).length : 1)
    const completeCreateNameParentLocalIdsByRenderedValues: CompleteCreateNameParentsByRenderedValues = useSelector(selectCompleteCreateNameParentLocalIdsByRenderedValues)
    const needsParentCompleteCreateName: boolean = useSelector(selectRuleParentLocalIdsByRuleName)[ruleName] != undefined
    let fullRenderedTokensPrefix: string | undefined = useSelector(selectRenderedCompleteCreateNamesByLocalId)[parentCompleteCreateNameLocalId]

    if (
        renderedTokensPrefix != undefined
        && parentCompleteCreateNameLocalId != undefined
    ) {
        fullRenderedTokensPrefix += renderedTokensPrefix

    } else if (!needsParentCompleteCreateName) {
        fullRenderedTokensPrefix = renderedTokensPrefix
    }

    const hasValue = value != undefined
    let counterValueCollision = false  // TODO

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

        const localMaxValue = hierarchyValues.length > 0 ? Math.max(...hierarchyValues) : -1

        fetchNextCounterValueFromMagma(projectName, ruleName, fullRenderedTokensPrefix)
            .then(remoteNextValue => {
                handleSetCounterValue(Math.max(localMaxValue + 1, remoteNextValue))
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
                    // disabled={hasValue && !counterValueCollision}
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