import {SessionStatusResponse} from "../../api_types";

export const statusWithDownloads: SessionStatusResponse = {
    "session": {
        "project_name": "labors",
        "workflow_name": "test_workflow.cwl",
        "key": "mykey",
        "inputs": {"someIntWithoutDefault": 123, "pickANum/num": 300}
    },
    "status": [[{
        "status": "complete",
        "name": "firstAdd",
        "downloads": {"sum": "https://vulcan.test/api/labors/data/e3338af9c88021783eb5f8c8162f00c2dac036e2/sum"}
    }, {
        "status": "complete",
        "name": "pickANum",
        "downloads": {"num": "https://vulcan.test/api/labors/data/f1a0d6b936caaa8fa30f02e91f97fe0173481a5f/num"}
    }, {
        "status": "complete",
        "name": "finalStep",
        "downloads": {"sum": "https://vulcan.test/api/labors/data/3a93800756007cc116476382433ad136f2dfad06/sum"}
    }, {"status": "complete", "name": "aPlot", "downloads": null}]],
    "outputs": {
        "status": "complete",
        "downloads": {"the_result": "https://vulcan.test/api/labors/data/9152441dfb808be8269617667755550aeee85a4d/the_result"}
    }
};