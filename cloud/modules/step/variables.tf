
variable "sfn_state_machine_name" {
  type        = string
  description = "The name of the state machine"
}

variable "sfn_state_machine_type" {
  default     = "STANDARD"
  description = "The type of state machine"
  validation {
    condition     = can(regex("^(STANDARD|EXPRESS)$", var.sfn_state_machine_type))
    error_message = "Must be STANDARD or EXPRESS."
  }
}

variable "sfn_state_machine_json_file" {
  type        = string
  description = "The states language definition for the state machine"
}

variable "sfn_state_machine_role_arn" {
  type        = string
  description = "The ARN for the role used by the state machine"
}

