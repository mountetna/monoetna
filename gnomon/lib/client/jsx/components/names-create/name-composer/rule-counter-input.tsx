import React from "react";
import { makeStyles } from '@material-ui/core/styles';
import { useSelector, useDispatch } from 'react-redux'
import ButtonBase from "@material-ui/core/ButtonBase";
import InputBase from "@material-ui/core/InputBase";
import FormHelperText from '@material-ui/core/FormHelperText';
import FormControl from '@material-ui/core/FormControl';
import { selectCounterValuesWithRuleName } from "../../../selectors/names";



const useStyles = (inputWidthEm: number) => makeStyles((theme) => ({
    ruleCounterField: {
        display: "inline-flex",
        flexDirection: "column",
        alignItems: "center",
    },
    ruleCounterInput: {
        fontWeight: "bold",
        width: `${inputWidthEm}em`,
        minWidth: "1.5em",
        maxWidth: "4em",
    },
    unsetOrError: {
        color: "red"
    },
    autoIncrementButton: {},
}));


const RuleCounterField = ({ ruleName, value, handleSetCounterValue }:
    { ruleName: string, value?: number, handleSetCounterValue: (value?: number) => void }) => {

    const classes = useStyles(value ? String(value).length : 1)()
    const counterValuesByRuleName: Record<string, Record<string, number>> = useSelector(selectCounterValuesWithRuleName)
    const ruleNameCounterValueCounts = counterValuesByRuleName[ruleName]
    const hasValue = value != undefined
    const counterValueCollision = (
        hasValue && value in ruleNameCounterValueCounts && ruleNameCounterValueCounts[value] > 1
    )


    // TODO: proper collision detection based on token values
    // ie: VAL1—not just RULE1
    const handleClickAutoIncrement = () => {
        const ruleNameCounterValues = Object.keys(ruleNameCounterValueCounts).map((key) => Number(key))
        const maxValue = ruleNameCounterValues.length ? Math.max(...ruleNameCounterValues) : -1
        handleSetCounterValue(maxValue + 1)
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