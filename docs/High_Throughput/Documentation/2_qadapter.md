```
[section.qadapter]
nodes = 1
ntasks_per_node = 20
walltime = 1-00:00:00
queue = partition_name
account = account_name
```

## Sections

### `[section.qadapter]`

A subsection of the high throughput workflows that specifies the queuing parameters for each job.

#### `nodes`

`int`: The number of nodes requested (Default defined in the `my_qadapter.yaml` file)

#### `ntasks_per_node`

`int`: The number of tasks per node (Default defined in the `my_qadapter.yaml` file)

#### `walltime`

`str`: The requested wall time in `DD-HH:MM:SS` format (Default defined in the `my_qadapter.yaml` file)

#### `queue`

`str`: The partition to submit the job to (Default defined in the `my_qadapter.yaml` file)

#### `account`

`str`: The account to charge for the calculation (Default defined in the `my_qadapter.yaml` file)
