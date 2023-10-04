import React from "react";
import { useSelector, useDispatch } from 'react-redux'
import { makeStyles } from '@material-ui/core/styles';
import FormGroup from '@material-ui/core/FormGroup';
import ButtonBase from "@material-ui/core/ButtonBase";
import Checkbox from "@material-ui/core/Checkbox";
import FileCopyOutlinedIcon from '@material-ui/icons/FileCopyOutlined';
import DeleteOutlineOutlinedIcon from "@material-ui/icons/DeleteOutlineOutlined";

import { CreateName, CreateNameGroup, CreateNameTokenValue, Rule, RuleParent, RuleToken, Token, TokenValue } from "../../../models";
import { selectRulesByName, selectTokenValuesByLocalId, selectTokens, selectRuleParentLocalIdsByRuleName, selectRuleParentsByLocalId, selectRuleTokenLocalIdsWithRuleName, selectRuleTokensByLocalId } from "../../../selectors/rules";
import { selectCreateNamesByLocalId, selectCreateNameWithLocalId, selectCreateNameLocalIdsWithGroupId, selectCreateNameTokenValueLocalIdsWithCreateNameLocalId, selectCreateNameTokenValuesByLocalId, selectCreateNameTokenValueLocalIdsByCreateNameLocalId, selectCreateNameLocalIdsByGroupId } from "../../../selectors/names";
import { addOrReplaceCreateNameTokenValue, setCreateNameRuleCounterValue, setCreateNameGroupsSelected, duplicateCreateNameGroup, deleteGroupsWithNames } from "../../../actions/names";
import { TokenSelect } from "./select";
import RuleCounterField from "./rule-counter-input";
import { createLocalId } from "../../../utils/models";


const useStyles = makeStyles((theme) => ({
    createNameElementsEditor: {
        display: "inline-flex"
    }
}));


const CreateNameElementsEditor = ({ createName, rule }: { createName: CreateName, rule: Rule }) => {
    const classes = useStyles()
    const dispatch = useDispatch()

    const tokens: Record<string, Token> = useSelector(selectTokens)

    const ruleTokensByLocalId: Record<string, RuleToken> = useSelector(selectRuleTokensByLocalId)
    const ruleTokens: RuleToken[] = useSelector(selectRuleTokenLocalIdsWithRuleName(rule.name))
        .map(rtLocalId => ruleTokensByLocalId[rtLocalId])

    const createNameTokenValuesByLocalId: Record<string, CreateNameTokenValue> = useSelector(selectCreateNameTokenValuesByLocalId)
    const createNameTokenValues: CreateNameTokenValue[] = (useSelector(selectCreateNameTokenValueLocalIdsWithCreateNameLocalId(createName.localId)) || [])
        .map(cntvLocalId => createNameTokenValuesByLocalId[cntvLocalId])

    const tokenValuesByLocalId: Record<string, TokenValue> = useSelector(selectTokenValuesByLocalId)

    // move this to <TokenSelect>?
    const setTokenValue = (value: TokenValue, ruleToken: RuleToken) => {
        const createNameTokenValue: CreateNameTokenValue = {
            localId: createLocalId(),
            tokenValueLocalId: value.name,
            createNameLocalId: createName.localId,
            ruleTokenLocalId: ruleToken.localId,
        }
        const oldCreateNameTokenValue = createNameTokenValues.find(
            cntv => cntv.ruleTokenLocalId == ruleToken.localId
        )
        dispatch(addOrReplaceCreateNameTokenValue(createNameTokenValue, oldCreateNameTokenValue))
    }

    const setRuleCounterValue = (value?: number) => {
        dispatch(setCreateNameRuleCounterValue(createName.localId, value))
    }

    return (
        // <span className="create-name-elements-editor">
        <FormGroup row className={classes.createNameElementsEditor}>
            {
                // copying array bc sort modifies in place
                [...ruleTokens].sort(rt => rt.ord).map(ruleToken => {
                    const createNameTokenValue = createNameTokenValues.find(
                        cntv => cntv.ruleTokenLocalId == ruleToken.localId
                    )
                    const currentValue = createNameTokenValue ?
                        tokenValuesByLocalId[createNameTokenValue.tokenValueLocalId] : undefined

                    return (
                        <React.Fragment key={ruleToken.localId}>
                            <TokenSelect
                                token={tokens[ruleToken.tokenName]}
                                value={currentValue}
                                onSetTokenValue={(value: TokenValue) => setTokenValue(value, ruleToken)}
                            ></TokenSelect>
                        </React.Fragment>
                    )
                })
            }
            {
                rule.hasCounter &&
                <RuleCounterField
                    ruleName={rule.name}
                    value={createName.ruleCounterValue}
                    handleSetCounterValue={setRuleCounterValue}
                />
            }
        </FormGroup>
    )
}


const CreateNameGroupComposer = ({ createNameGroup }: { createNameGroup: CreateNameGroup }) => {
    const dispatch = useDispatch()
    
    const createNameLocalIds: string[] = useSelector(selectCreateNameLocalIdsWithGroupId(createNameGroup.localId))
    const createNamesByLocalId: Record<string, CreateName> = useSelector(selectCreateNamesByLocalId)
    const primaryCreateName: CreateName = useSelector(selectCreateNameWithLocalId(createNameGroup.primaryCreateNameId))
    const createNameLocalIdsByCreateNameGroupLocalId: Record<string, string[]> = useSelector(selectCreateNameLocalIdsByGroupId)
    const createNameTokenValueLocalIdsByCreateNameLocalId: Record<string, string[]> = useSelector(selectCreateNameTokenValueLocalIdsByCreateNameLocalId)
    const createNameTokenValuesByLocalId: Record<string, CreateNameTokenValue> = useSelector(selectCreateNameTokenValuesByLocalId)
    const allRules: Record<string, Rule> = useSelector(selectRulesByName)
    const ruleParentLocalIdsbyRuleName: Record<string, string[]> = useSelector(selectRuleParentLocalIdsByRuleName)
    const ruleParentsbyLocalId: Record<string, RuleParent> = useSelector(selectRuleParentsByLocalId)

    const orderedRuleNames: string[] = []

    const ruleNamesToScan: string[] = [primaryCreateName.ruleName]
    while (ruleNamesToScan.length) {
        const ruleName = ruleNamesToScan.pop()
        if (!ruleName) {
            break
        }
        const rule = allRules[ruleName]

        orderedRuleNames.unshift(rule.name)

        if (ruleParentLocalIdsbyRuleName[ruleName]) {
            const parentRuleNames = ruleParentLocalIdsbyRuleName[ruleName].map(rpLocalId => {
                return ruleParentsbyLocalId[rpLocalId].parentRuleName
            })
            ruleNamesToScan.push(...parentRuleNames)
        }
    }

    const handleClickSelect = (event: React.ChangeEvent) => {
        dispatch(setCreateNameGroupsSelected([createNameGroup.localId], event.target.checked))
    }

    const handleClickCopy = () => {
        dispatch(duplicateCreateNameGroup(
            createNameGroup,
            createNameLocalIdsByCreateNameGroupLocalId,
            createNamesByLocalId,
            createNameTokenValueLocalIdsByCreateNameLocalId,
            createNameTokenValuesByLocalId,
        ))
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
                    const createNameLocalId = createNameLocalIds.find((cnLocalId) => {
                        return createNamesByLocalId[cnLocalId].ruleName == ruleName
                    })
                    if (!createNameLocalId) {
                        console.error(`Error creating CreateNameElementsEditor. CreateName with ruleName ${ruleName} not found.`)
                        return
                    }
                    const rule = allRules[ruleName]

                    return (
                        <React.Fragment key={rule.name}>
                            <CreateNameElementsEditor
                                createName={createNamesByLocalId[createNameLocalId]}
                                rule={rule}
                            />
                        </React.Fragment>
                    )
                })
            }
        </div>
    )
}


export default CreateNameGroupComposer;