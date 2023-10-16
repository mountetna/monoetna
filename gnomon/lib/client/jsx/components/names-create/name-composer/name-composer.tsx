import React from "react";
import { useSelector, useDispatch } from 'react-redux'
import { makeStyles } from '@material-ui/core/styles';
import FormGroup from '@material-ui/core/FormGroup';
import ButtonBase from "@material-ui/core/ButtonBase";
import Checkbox from "@material-ui/core/Checkbox";
import FileCopyOutlinedIcon from '@material-ui/icons/FileCopyOutlined';
import DeleteOutlineOutlinedIcon from "@material-ui/icons/DeleteOutlineOutlined";

import { CreateName, CreateNameGroup, CreateNameTokenValue, Rule, RuleParent, RuleToken, Token, TokenValue, UNSET_TOKEN_VALUE, UNSET_VALUE } from "../../../models";
import { selectRulesByName, selectTokenValuesByLocalId, selectTokens, selectRuleParentLocalIdsByRuleName, selectRuleParentsByLocalId, selectRuleTokenLocalIdsWithRuleName, selectRuleTokensByLocalId, selectRulesNamesHierarchicalListByPrimaryRuleName } from "../../../selectors/rules";
import { selectCreateNamesByLocalId, selectCreateNameWithLocalId, selectCreateNameLocalIdsWithGroupId, selectCreateNameTokenValueLocalIdsWithCreateNameLocalId, selectCreateNameTokenValuesByLocalId, selectCreateNameTokenValueLocalIdsByCreateNameLocalId, selectCreateNameLocalIdsByGroupId, selectSelectedCreateNameGroupIds } from "../../../selectors/names";
import { addOrReplaceCreateNameTokenValues, setCreateNameRuleCounterValues, duplicateCreateNameGroups, deleteGroupsWithNames, addCreateNameGroupsToSelection, removeCreateNameGroupsFromSelection } from "../../../actions/names";
import { TokenSelect } from "./select";
import RuleCounterField from "./rule-counter-input";
import { createLocalId } from "../../../utils/models";


const useStyles = makeStyles((theme) => ({
    createNameElementsEditor: {
        display: "inline-flex"
    }
}));


const CreateNameElementsEditor = ({ createName, rule, includeUnsetAsValue }: { createName: CreateName, rule: Rule, includeUnsetAsValue: boolean }) => {
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

    const setTokenValue = (value: TokenValue, ruleToken: RuleToken) => {
        const createNameTokenValue: CreateNameTokenValue = {
            localId: createLocalId(),
            tokenValueLocalId: value.localId,
            createNameLocalId: createName.localId,
            ruleTokenLocalId: ruleToken.localId,
        }

        if (value == UNSET_TOKEN_VALUE) {
            createNameTokenValue.tokenValueLocalId = UNSET_VALUE
        }

        const oldCreateNameTokenValue = createNameTokenValues.find(
            cntv => cntv.ruleTokenLocalId == ruleToken.localId
        )
        dispatch(addOrReplaceCreateNameTokenValues([createNameTokenValue], oldCreateNameTokenValue ? [oldCreateNameTokenValue.localId] : []))
    }

    const setRuleCounterValue = (value?: number) => {
        dispatch(setCreateNameRuleCounterValues({ [createName.localId]: value }))
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
                                includeUnsetAsValue={includeUnsetAsValue}
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


const CreateNameGroupComposer = ({ createNameGroup, includeTools = false, includeUnsetAsValue = false }: {
    createNameGroup: CreateNameGroup,
    includeTools?: boolean
    includeUnsetAsValue?: boolean
}) => {
    const dispatch = useDispatch()

    const createNameLocalIds: string[] = useSelector(selectCreateNameLocalIdsWithGroupId(createNameGroup.localId))
    const createNamesByLocalId: Record<string, CreateName> = useSelector(selectCreateNamesByLocalId)
    const primaryCreateName: CreateName = useSelector(selectCreateNameWithLocalId(createNameGroup.primaryCreateNameLocalId))
    const createNameLocalIdsByCreateNameGroupLocalId: Record<string, string[]> = useSelector(selectCreateNameLocalIdsByGroupId)
    const createNameTokenValueLocalIdsByCreateNameLocalId: Record<string, string[]> = useSelector(selectCreateNameTokenValueLocalIdsByCreateNameLocalId)
    const createNameTokenValuesByLocalId: Record<string, CreateNameTokenValue> = useSelector(selectCreateNameTokenValuesByLocalId)
    const allRules: Record<string, Rule> = useSelector(selectRulesByName)
    const orderedRuleNames: string[] = useSelector(selectRulesNamesHierarchicalListByPrimaryRuleName)[primaryCreateName.ruleName]
    const selectedCreateNameGroupsIds: Set<string> = useSelector(selectSelectedCreateNameGroupIds)

    const handleClickSelect = (event: React.ChangeEvent) => {
        if (event.target.checked) {
            dispatch(addCreateNameGroupsToSelection([createNameGroup.localId]))
            return
        }
        dispatch(removeCreateNameGroupsFromSelection([createNameGroup.localId]))
    }

    const handleClickCopy = () => {
        dispatch(duplicateCreateNameGroups(
            [createNameGroup],
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
            {
                includeTools &&
                <span className="create-name-group-composer-tools">
                    <Checkbox
                        checked={selectedCreateNameGroupsIds.has(createNameGroup.localId)}
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
            }
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
                                includeUnsetAsValue={includeUnsetAsValue}
                            />
                        </React.Fragment>
                    )
                })
            }
        </div>
    )
}


export default CreateNameGroupComposer;