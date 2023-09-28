import React from "react";
import { useSelector, useDispatch } from 'react-redux'
import { makeStyles } from '@material-ui/core/styles';
import FormGroup from '@material-ui/core/FormGroup';
import ButtonBase from "@material-ui/core/ButtonBase";
import Checkbox from "@material-ui/core/Checkbox";
import FileCopyOutlinedIcon from '@material-ui/icons/FileCopyOutlined';
import DeleteOutlineOutlinedIcon from "@material-ui/icons/DeleteOutlineOutlined";

import { CreateName, CreateNameGroup, Rule, TOKEN_VALUE_PLACEHOLDER, Token, TokenValue } from "../../../models";
import { selectRules, selectTokens } from "../../../selectors/rules";
import { selectCreateNameById, selectCreateNames, selectCreateNamesByIds } from "../../../selectors/names";
import { setCreateNameTokenValue, setCreateNameCounterValue, setCreateNameGroupsSelected, duplicateCreateNameGroup, deleteGroupsWithNames } from "../../../actions/names";
import TokenSelect from "./token-select";
import RuleCounterField from "./rule-counter-input";


const useStyles = makeStyles((theme) => ({
    createNameElementsEditor: {
        display: "inline-flex"
    }
}));


const CreateNameElementsEditor = ({ createName, rule }: { createName: CreateName, rule: Rule }) => {
    const classes = useStyles()
    const tokens: Record<string, Token> = useSelector(selectTokens)
    const dispatch = useDispatch()

    const setTokenValue = (value: TokenValue, idx: number) => {
        dispatch(setCreateNameTokenValue(createName.localId, value, idx))
    }

    const setRuleCounterValue = (value?: number) => {
        dispatch(setCreateNameCounterValue(createName.localId, value))
    }

    return (
        // <span className="create-name-elements-editor">
        <FormGroup row className={classes.createNameElementsEditor}>
            {
                rule.tokenNames.map((tokenName, idx) => {
                    const currentValue = createName.tokenValues[idx]

                    return (
                        <React.Fragment key={idx}>
                            <TokenSelect
                                token={tokens[tokenName]}
                                value={currentValue != TOKEN_VALUE_PLACEHOLDER ? currentValue : undefined}
                                handleSetTokenValue={(value: TokenValue) => setTokenValue(value, idx)}
                            ></TokenSelect>
                        </React.Fragment>
                    )
                })
            }
            {
                rule.hasCounter &&
                <RuleCounterField
                    ruleName={rule.name}
                    value={createName.counterValue}
                    handleSetCounterValue={setRuleCounterValue}
                />
            }
        </FormGroup>
    )
}


const CreateNameGroupComposer = ({ createNameGroup }: { createNameGroup: CreateNameGroup }) => {
    const dispatch = useDispatch()
    const createNamesById: Record<string, CreateName> = useSelector(selectCreateNames)
    const createNames: CreateName[] = useSelector(selectCreateNamesByIds(createNameGroup.createNameIds))
    const primaryCreateName: CreateName = useSelector(selectCreateNameById(createNameGroup.primaryCreateNameId))
    const allRules: Record<string, Rule> = useSelector(selectRules)
    const orderedRuleNames: string[] = []

    const ruleNamesToScan: string[] = [primaryCreateName.ruleName]
    while (ruleNamesToScan.length) {
        const ruleName = ruleNamesToScan.pop()
        if (!ruleName) {
            break
        }
        const rule = allRules[ruleName]

        orderedRuleNames.unshift(rule.name)

        if (rule.parentRuleNames.length) {
            ruleNamesToScan.push(...rule.parentRuleNames)
        }
    }

    const handleClickSelect = (event: React.ChangeEvent) => {
        dispatch(setCreateNameGroupsSelected([createNameGroup.localId], event.target.checked))
    }

    const handleClickCopy = () => {
        dispatch(duplicateCreateNameGroup(createNameGroup, createNamesById))
    }

    const handleClickDelete = () => {
        dispatch(deleteGroupsWithNames([createNameGroup.localId]))
    }

    return (
        <div className="create-name-group-composer">
            <span className="create-name-group-composer-tools">
                <Checkbox
                    checked={createNameGroup.selected}
                    onChange={handleClickSelect}
                    inputProps={{ 'aria-label': 'Select Name' }}
                />
                <ButtonBase
                    onClick={handleClickCopy}
                    aria-label={`Copy Name with Values"`}
                    disableRipple
                    disableTouchRipple
                >
                    <FileCopyOutlinedIcon />
                </ButtonBase>
                <ButtonBase
                    onClick={handleClickDelete}
                    aria-label={`Delete Name"`}
                    disableRipple
                    disableTouchRipple
                >
                    <DeleteOutlineOutlinedIcon />
                </ButtonBase>
            </span>
            {
                orderedRuleNames.map((ruleName) => {
                    const createName = createNames.find((createName) => createName.ruleName == ruleName)
                    if (!createName) {
                        console.error(`Error creating CreateNameElementsEditor. CreateName with ruleName ${ruleName} not found.`)
                        return
                    }
                    const rule = allRules[ruleName]

                    return (
                        <React.Fragment key={rule.name}>
                            <CreateNameElementsEditor createName={createName} rule={rule} />
                        </React.Fragment>
                    )
                })
            }
        </div>
    )
}


export default CreateNameGroupComposer;