import React, { useEffect, useState } from "react";
import { useSelector } from 'react-redux'
import { makeStyles } from '@material-ui/core/styles';
import FormGroup from '@material-ui/core/FormGroup';
import ButtonBase from "@material-ui/core/ButtonBase";
import Checkbox from "@material-ui/core/Checkbox";
import FileCopyOutlinedIcon from '@material-ui/icons/FileCopyOutlined';
import DeleteOutlineOutlinedIcon from "@material-ui/icons/DeleteOutlineOutlined";
import _ from "lodash"

import { CreateName, CreateNameGroup, CreateNameTokenValue, Rule, RuleToken, TokenValue, UNSET_TOKEN_VALUE, UNSET_VALUE } from "../../../models";
import { selectRulesByName, selectTokenValuesByLocalId, selectTokens, selectRuleTokenLocalIdsWithRuleName, selectRuleTokensByLocalId, selectRuleNamesHierarchicalListByPrimaryRuleName } from "../../../selectors/rules";
import { selectCreateNamesByLocalId, selectCreateNameWithLocalId, selectCreateNameLocalIdsWithGroupId, selectCreateNameTokenValueLocalIdsWithCreateNameLocalId, selectCreateNameTokenValuesByLocalId, selectSelectedCreateNameGroupIds, selectSortedCompleteCreateNamesWithCreateNameGroupLocalId, selectCreateNameCompleteCreateNameLocalIdsByCompleteCreateNameLocalId, selectRenderedCompleteCreateNamesByLocalId } from "../../../selectors/names";
import { addOrReplaceCreateNameTokenValues, setCreateNameRuleCounterValues, duplicateCreateNameGroups, deleteGroupsWithNames, addCreateNameGroupsToSelection, removeCreateNameGroupsFromSelection, deleteCreateNameTokenValue } from "../../../actions/names";
import { selectPathParts } from "../../../selectors/location"
import { TokenSelect } from "./select";
import RuleCounterField from "./rule-counter-input";
import { createLocalId } from "../../../utils/models";
import { useDispatch } from "../../../utils/redux";
import { selectGlobalState } from "../../../selectors/global";
import { fetchWhetherNameExistsInMagma } from "../../../utils/names";


const useEditorStyles = makeStyles((theme) => ({
    elementsEditorContainer: {
        display: "inline-flex",
        alignItems: "center",
        flexWrap: "nowrap",
        // account for absolute positioning of <RuleCounterField>
        marginRight: "0.4em",
    },
}));


const CreateNameElementsEditor = ({ createName, rule, includeUnsetAsValue, parentCompleteCreateNameLocalId, error }: {
    createName: CreateName,
    rule: Rule,
    includeUnsetAsValue: boolean,
    parentCompleteCreateNameLocalId: string | undefined,
    error?: boolean
}) => {
    const classes = useEditorStyles()
    const dispatch = useDispatch()

    const globalState = useSelector(selectGlobalState)
    const projectName = useSelector(selectPathParts)[0]
    const tokens = useSelector(selectTokens)

    const ruleTokensByLocalId = useSelector(selectRuleTokensByLocalId)
    let sortedRuleTokens = useSelector(selectRuleTokenLocalIdsWithRuleName(rule.name))
        .map(rtLocalId => ruleTokensByLocalId[rtLocalId])

    // copying array bc sort modifies in place
    sortedRuleTokens = [...sortedRuleTokens].sort(rt => rt.ord)

    const createNameTokenValuesByLocalId = useSelector(selectCreateNameTokenValuesByLocalId)
    const createNameTokenValues = (useSelector(selectCreateNameTokenValueLocalIdsWithCreateNameLocalId(createName.localId)) || [])
        .map(cntvLocalId => createNameTokenValuesByLocalId[cntvLocalId])

    const tokenValuesByLocalId = useSelector(selectTokenValuesByLocalId)

    const sortedTokenValues = sortedRuleTokens.map(ruleToken => {
        const createNameTokenValue = createNameTokenValues.find(
            cntv => cntv.ruleTokenLocalId == ruleToken.localId
        )
        return createNameTokenValue ?
            tokenValuesByLocalId[createNameTokenValue.tokenValueLocalId] : undefined
    })

    const renderedTokens: string | undefined = _.every(sortedTokenValues, tv => tv != undefined)
        ? sortedTokenValues.map(tv => tv.name).reduce((prev, curr) => prev + curr, "")
        : undefined


    const setTokenValue = (ruleToken: RuleToken, value?: TokenValue) => {
        const oldCreateNameTokenValue = createNameTokenValues.find(
            cntv => cntv.ruleTokenLocalId == ruleToken.localId
        )

        if (!value) {
            if (oldCreateNameTokenValue) {
                dispatch(deleteCreateNameTokenValue(oldCreateNameTokenValue, globalState))
            }
            return
        }

        const createNameTokenValue: CreateNameTokenValue = {
            localId: createLocalId(),
            tokenValueLocalId: value.localId,
            createNameLocalId: createName.localId,
            ruleTokenLocalId: ruleToken.localId,
        }

        if (value == UNSET_TOKEN_VALUE) {
            createNameTokenValue.tokenValueLocalId = UNSET_VALUE
        }

        dispatch(
            addOrReplaceCreateNameTokenValues(
                [createNameTokenValue],
                oldCreateNameTokenValue ? [oldCreateNameTokenValue.localId] : [],
                globalState,
            )
        )
    }

    const setRuleCounterValue = (value?: number) => {
        dispatch(
            setCreateNameRuleCounterValues(
                { [createName.localId]: value },
                globalState,
            )
        )
    }

    return (
        <FormGroup row className={classes.elementsEditorContainer}>
            {
                sortedRuleTokens.map((ruleToken, idx) =>

                    <React.Fragment key={ruleToken.localId}>
                        <TokenSelect
                            token={tokens[ruleToken.tokenName]}
                            value={sortedTokenValues[idx]}
                            onSetTokenValue={value => setTokenValue(ruleToken, value)}
                            includeUnsetAsValue={includeUnsetAsValue}
                        ></TokenSelect>
                    </React.Fragment>
                )
            }
            {
                rule.hasCounter &&
                <RuleCounterField
                    value={createName.ruleCounterValue}
                    renderedTokensPrefix={renderedTokens}
                    parentCompleteCreateNameLocalId={parentCompleteCreateNameLocalId}
                    projectName={projectName}
                    ruleName={rule.name}
                    handleSetCounterValue={setRuleCounterValue}
                />
            }
        </FormGroup>
    )
}


const useComposerStyles = makeStyles((theme) => ({
    container: {
        display: "inline-flex",
        flexWrap: "wrap",
        "& > *": {
            // account for absolute positioning of <RuleCounterField>
            paddingTop: "1.5em",
        },
        "&.hasError > *": {
            padding: "1.5em 0",
        }
    },
    toolsContainer: {
        display: "inline-flex",
        alignItems: "center",
        marginRight: "1em",

    },
    checkbox: {
        padding: "0",
        paddingRight: "0.25em"
    },
    editorsContainer: {},
    errorMessage: {
        position: "absolute",
        color: "red",
        fontStyle: "italic",
        fontSize: "14px",
    },
}));


const CreateNameGroupComposer = ({ createNameGroup, className, includeTools = false, includeUnsetAsValue = false }: {
    createNameGroup: CreateNameGroup,
    className?: string,
    includeTools?: boolean
    includeUnsetAsValue?: boolean
}) => {
    const dispatch = useDispatch()
    const classes = useComposerStyles()
    const [duplicateTracker, setDuplicateTracker] = useState({ local: 0, remote: 0 })

    const globalState = useSelector(selectGlobalState)
    const projectName = useSelector(selectPathParts)[0]
    const createNameLocalIds = useSelector(selectCreateNameLocalIdsWithGroupId(createNameGroup.localId))
    const createNamesByLocalId = useSelector(selectCreateNamesByLocalId)
    const primaryCreateName = useSelector(selectCreateNameWithLocalId(createNameGroup.primaryCreateNameLocalId))
    const allRules = useSelector(selectRulesByName)
    const orderedRuleNames = useSelector(selectRuleNamesHierarchicalListByPrimaryRuleName)[primaryCreateName.ruleName]
    const selectedCreateNameGroupsIds = useSelector(selectSelectedCreateNameGroupIds)

    const sortedCompleteCreateNamesWithCreateNameGroupLocalId = useSelector(
        selectSortedCompleteCreateNamesWithCreateNameGroupLocalId(createNameGroup.localId)
    )
    const primaryCompleteCreateName = sortedCompleteCreateNamesWithCreateNameGroupLocalId[sortedCompleteCreateNamesWithCreateNameGroupLocalId.length - 1]
    const createNameCompleteCreateNameLocalIdsByCompleteCreateNameLocalId = useSelector(selectCreateNameCompleteCreateNameLocalIdsByCompleteCreateNameLocalId)
    const renderedCompleteCreateNamesByLocalId = useSelector(selectRenderedCompleteCreateNamesByLocalId)

    useEffect(() => {
        if (!primaryCompleteCreateName) {
            setDuplicateTracker({ local: 0, remote: 0 })
            return
        }

        const localDuplicates = createNameCompleteCreateNameLocalIdsByCompleteCreateNameLocalId[
            primaryCompleteCreateName.localId
        ].length

        const renderedName = renderedCompleteCreateNamesByLocalId[primaryCompleteCreateName.localId]

        fetchWhetherNameExistsInMagma(projectName, primaryCreateName.ruleName, renderedName)
            .then(nameExists => {
                setDuplicateTracker({ local: localDuplicates, remote: nameExists ? 1 : 0 })
            })
            .catch(err => `Error determining whether name "${renderedName}" has remote duplicate: ${err}"`)

    }, [primaryCompleteCreateName?.localId])

    const createErrorMessage = () => {
        const errorMsgs: string[] = []

        if (duplicateTracker.local > 1) {
            errorMsgs.push("locally")
        }
        if (duplicateTracker.remote) {
            errorMsgs.push("in Magma")
        }

        return errorMsgs.length
            ? `Name already exists ${errorMsgs.join(" and ")}`
            : undefined
    }

    const errorMessage = createErrorMessage()

    const handleClickSelect = (event: React.ChangeEvent<HTMLInputElement>) => {
        if (event.target.checked) {
            dispatch(addCreateNameGroupsToSelection([createNameGroup.localId]))
            return
        }
        dispatch(removeCreateNameGroupsFromSelection([createNameGroup.localId]))
    }

    const handleClickCopy = () => {
        dispatch(duplicateCreateNameGroups([createNameGroup], globalState))
    }

    const handleClickDelete = () => {
        dispatch(deleteGroupsWithNames([createNameGroup.localId], globalState))
    }

    return (
        <div
            className={`${classes.container} ${errorMessage ? "hasError" : ""} ${className != undefined ? className : ""}`}
        >
            {
                includeTools &&
                <span className={classes.toolsContainer}>
                    <Checkbox
                        checked={selectedCreateNameGroupsIds.has(createNameGroup.localId)}
                        onChange={handleClickSelect}
                        inputProps={{ 'aria-label': 'Select Name' }}
                        className={classes.checkbox}
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
            <span className={classes.editorsContainer}>
                {
                    orderedRuleNames.map((ruleName, idx) => {
                        const createNameLocalId = createNameLocalIds.find((cnLocalId) => {
                            return createNamesByLocalId[cnLocalId].ruleName == ruleName
                        })
                        if (!createNameLocalId) {
                            console.error(`Error creating CreateNameElementsEditor. CreateName with ruleName ${ruleName} not found.`)
                            return
                        }
                        const rule = allRules[ruleName]
                        const parentCompleteCreateName = sortedCompleteCreateNamesWithCreateNameGroupLocalId[idx - 1]

                        return (
                            <React.Fragment key={rule.name}>
                                <CreateNameElementsEditor
                                    createName={createNamesByLocalId[createNameLocalId]}
                                    rule={rule}
                                    parentCompleteCreateNameLocalId={parentCompleteCreateName?.localId}
                                    includeUnsetAsValue={includeUnsetAsValue}
                                    error={errorMessage != undefined}
                                />
                            </React.Fragment>
                        )
                    })
                }
                {errorMessage &&
                    <div className={classes.errorMessage}>{errorMessage}</div>}
            </span>
        </div>
    )
}


export default CreateNameGroupComposer;