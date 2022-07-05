import pytest

@pytest.fixture
def labor_template():
    return {
        "name": "labor",
        "attributes": {
            "created_at": {
                "name": "created_at",
                "attribute_name": "created_at",
                "display_name": "Created At",
                "hidden": True,
                "validation": None,
                "attribute_type": "date_time"
            },
            "updated_at": {
                "name": "updated_at",
                "attribute_name": "updated_at",
                "display_name": "Updated At",
                "hidden": True,
                "validation": None,
                "attribute_type": "date_time"
            },
            "name": {
                "name": "name",
                "attribute_name": "name",
                "description": "Name for this labor.",
                "display_name": "Name",
                "restricted": False,
                "read_only": False,
                "hidden": False,
                "validation": None,
                "attribute_type": "identifier"
            },
            "country": {
                "name": "country",
                "attribute_name": "country",
                "description": "Place where this labor took place.",
                "restricted": False,
                "read_only": False,
                "hidden": False,
                "validation": None,
                "attribute_type": "string"
            }
        },
        "identifier": "name"
    }
